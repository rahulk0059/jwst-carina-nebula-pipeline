"""Step 2 of the per-filter procedure: measure the raw gauge against an official i2d.

The cross-visit correction applied by ``mosaic`` is a *relative* translation, so
its absolute origin is a free choice and nothing in the star matches identifies
it.  This module measures that origin: it star-detects one reference frame per
``(visit, detector)`` group, converts those stars into the pixel frame of the
official combined i2d, and reports the median ``dx``/``dy`` of each group versus
the i2d.  The visit whose offset is nearest zero is the anchor, and that is the
value to pass as ``mosaic --gauge-visit``.

The offsets reported here are the **raw, uncorrected** header-WCS error, so
they are a diagnostic and not an input to the solve.  Gauging the mosaic to the
i2d instead of to a measured anchor would drive the validated residual to ~0 by
construction and destroy the evidence; see "The gauge is the part that is easy
to get wrong" in AGENTS.md.

Units are i2d pixels throughout, and the match radius is specified in arcsec so
that it means the same angular distance in every channel.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from astropy.io import fits

from jwst_stack.io import CalExposure, mask_scaled_sci
from jwst_stack.register import detect_stars, match_stars
from jwst_stack.units import I2D_PX, px_to_arcsec
from jwst_stack.validation import open_i2d

#: NIRCam samples its PSF at roughly 2-3 px FWHM in *every* channel, so this is
#: filter-independent even though the pixel scale is not.  ``stack``/``mosaic``
#: default to 3.0; this module uses 2.5, which is what the recorded F187N and
#: F335M step-2 runs used, so re-running them reproduces their numbers.
DEFAULT_FWHM_PX = 2.5

#: Must exceed the largest offset being measured (inter-visit pointing runs to
#: ~1 px) while staying far below the ~32 px at short-wave scale that
#: cross-matches unrelated stars in this crowded field.  0.3 arcsec is 4.8 px
#: in the long-wave channel and 9.7 px in the short-wave one.
DEFAULT_MATCH_ARCSEC = 0.3

#: A group with fewer matches than this has no usable median; it is reported as
#: skipped rather than folded into the per-visit median.
DEFAULT_MIN_MATCHED = 20

#: F335M and F187N both cleared this floor with 109-604 matches per group.
DEFAULT_MAX_I2D_STARS = 8000
DEFAULT_MAX_FRAME_STARS = 400

#: A group's median match distance above this fraction of the match radius
#: means the offset is being measured at the edge of the search window, where a
#: real translation and a lucky pairing with a neighbouring star are not
#: distinguishable.  Every shipped filter sits near 10% of its radius
#: (F335M 0.49 px of 4.77, F187N/F200W/F090W 0.9 px of 9.7), so this only
#: fires when the radius is too small for the offset being measured.
#:
#: It deliberately does **not** try to detect "wrong" matches via scatter.  A
#: regular star lattice shifted past the radius pairs each star with a
#: neighbour and the resulting median is tight to ~0.02 px, so a scatter test
#: passes a measurement that is entirely spurious.  There is no cheap
#: statistic that separates the two cases: both produce the same median match
#: distance.  The diagnostic that does help is the match *rate*, which the
#: report prints next to the radius.
SUSPECT_RADIUS_FRACTION = 0.4


class GaugeError(RuntimeError):
    """Raised when the gauge cannot be measured at all (as opposed to skipped groups)."""


@dataclass
class GaugeGroup:
    """Raw offset of one ``(visit, detector)`` group versus the official i2d."""

    visit: int
    detector: str
    n_matched: int
    dx_median: float
    dy_median: float
    offset_median: float
    offset_mad: float
    #: Frame stars offered, before matching.  The ratio to ``n_matched`` is the
    #: diagnostic that separates a real offset from noise: a frame whose stars
    #: genuinely match the i2d loses a handful, while a frame that is matched to
    #: the wrong stars typically pairs only a fraction of them.  Recorded
    #: because it is the one number here that can expose a bad match, and the
    #: F090W/F200W step-2 JSONs already carry extra fields beyond F187N's seven.
    n_frame_stars: int = 0

    @property
    def match_rate(self) -> float:
        """Fraction of the frame's detected stars that matched, or nan if none."""
        if self.n_frame_stars <= 0:
            return float("nan")
        return self.n_matched / self.n_frame_stars


@dataclass
class GaugeVisit:
    """Per-visit roll-up: the median of that visit's per-detector medians."""

    visit: int
    dx_median: float
    dy_median: float
    offset_median: float
    spread: float
    n_det: int
    internal_rms: float


@dataclass
class GaugeResult:
    """Everything :func:`measure_gauge` learned, in i2d pixels."""

    groups: dict[str, GaugeGroup] = field(default_factory=dict)
    visits: dict[int, GaugeVisit] = field(default_factory=dict)
    anchor_visit: int | None = None
    i2d_scale_arcsec: float = 0.0
    i2d_shape: tuple[int, int] = (0, 0)
    match_arcsec: float = DEFAULT_MATCH_ARCSEC
    n_i2d_stars: int = 0
    skipped: list[str] = field(default_factory=list)

    @property
    def offset_median(self) -> float:
        """Offset of the anchor visit, i.e. how far the gauge sits from the i2d."""
        if self.anchor_visit is None:
            return float("nan")
        return self.visits[self.anchor_visit].offset_median

    @property
    def offset_median_arcsec(self) -> float:
        return self.offset_median * self.i2d_scale_arcsec

    def separation_from_anchor(self, visit: int) -> float:
        """Visit-to-anchor separation in i2d px, the quantity to predict from."""
        if self.anchor_visit is None:
            return float("nan")
        a = self.visits[self.anchor_visit]
        v = self.visits[visit]
        return float(np.hypot(v.dx_median - a.dx_median, v.dy_median - a.dy_median))

    def as_groups_dict(self) -> dict[str, dict]:
        """Flat ``{"v1_nrca1": {...}}`` form, matching the shipped step-2 JSONs."""
        return {k: asdict(v) for k, v in self.groups.items()}

    def as_dict(self) -> dict:
        return {
            "groups": self.as_groups_dict(),
            "visits": {str(v): asdict(s) for v, s in self.visits.items()},
            "anchor_visit": self.anchor_visit,
            "i2d_scale_arcsec": self.i2d_scale_arcsec,
            "i2d_shape": list(self.i2d_shape),
            "match_arcsec": self.match_arcsec,
            "n_i2d_stars": self.n_i2d_stars,
            "skipped": list(self.skipped),
        }


def _i2d_stars(
    i2d_path: str | Path,
    fwhm_px: float,
    threshold_sigma: float,
    max_stars: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, tuple[int, int]]:
    """Detect the reference star list in the official i2d, in i2d pixels."""
    sci, wcs, shape = open_i2d(i2d_path)
    scale = float(wcs.proj_plane_pixel_scales()[0].to_value("arcsec"))
    image = mask_scaled_sci(np.asarray(sci, dtype=np.float32), None)
    with fits.open(i2d_path) as hdul:
        if "WHT" in hdul:
            # Uncovered mosaic area is exactly zero; without this the star
            # finder runs on empty background and the noise estimate is set by
            # it rather than by the covered pixels.
            image[np.asarray(hdul["WHT"].data) <= 0] = np.nan
    x, y, flux = detect_stars(
        image, fwhm_px=fwhm_px, threshold_sigma=threshold_sigma, max_stars=max_stars
    )
    return x, y, flux, scale, shape


def measure_gauge(
    exposures: Iterable[CalExposure],
    i2d_path: str | Path,
    *,
    reference_exposure: int = 1,
    match_arcsec: float = DEFAULT_MATCH_ARCSEC,
    fwhm_px: float = DEFAULT_FWHM_PX,
    threshold_sigma: float = 5.0,
    max_i2d_stars: int = DEFAULT_MAX_I2D_STARS,
    max_frame_stars: int = DEFAULT_MAX_FRAME_STARS,
    min_matched: int = DEFAULT_MIN_MATCHED,
) -> GaugeResult:
    """Measure every ``(visit, detector)`` group's raw offset from the i2d.

    One frame per group -- the lowest ``exposure`` number, or *reference_exposure*
    if given -- supplies the stars, matching the recorded step-2 runs.  Groups
    with fewer than *min_matched* matches are listed in
    :attr:`GaugeResult.skipped` rather than silently averaged in; a skipped group
    is not an error, but it does weaken that visit's median, and the visit's
    ``n_det`` says so.

    Raises :class:`GaugeError` only when nothing can be measured, since a gauge
    that cannot be read must not be mistaken for a gauge of zero.
    """
    groups: dict[str, GaugeGroup] = {}
    skipped: list[str] = []

    ix, iy, _, scale, (ny, nx) = _i2d_stars(
        i2d_path, fwhm_px, threshold_sigma, max_i2d_stars
    )
    if len(ix) < min_matched:
        raise GaugeError(
            f"only {len(ix)} stars detected in {i2d_path}; need >= {min_matched}"
        )

    _, i2d_wcs, _ = open_i2d(i2d_path)
    radius_px = match_arcsec / scale

    for exp in exposures:
        if exp.exposure != reference_exposure:
            continue
        key = f"v{exp.visit}_{exp.detector.strip().lower()}"
        if key in groups:
            continue

        frame = mask_scaled_sci(exp.sci, exp.dq)
        fx, fy, _ = detect_stars(
            frame, fwhm_px=fwhm_px, threshold_sigma=threshold_sigma,
            max_stars=max_frame_stars,
        )
        if len(fx) < min_matched:
            skipped.append(f"{key}: {len(fx)} frame stars")
            continue

        ra, dec = exp.wcs.all_pix2world(fx, fy, 0)
        gx, gy = i2d_wcs.all_world2pix(np.atleast_1d(ra), np.atleast_1d(dec), 0)
        gx = np.asarray(gx, dtype=float)
        gy = np.asarray(gy, dtype=float)
        inside = (
            np.isfinite(gx) & np.isfinite(gy)
            & (gx >= 0) & (gx < nx) & (gy >= 0) & (gy < ny)
        )
        dx, dy, resid = match_stars(gx[inside], gy[inside], ix, iy, radius_px)

        if len(dx) < min_matched:
            skipped.append(f"{key}: {len(dx)} matched of {len(fx)} frame stars")
            continue

        offset = float(np.median(resid))
        groups[key] = GaugeGroup(
            visit=exp.visit,
            detector=exp.detector.strip().lower(),
            n_matched=int(len(dx)),
            dx_median=float(np.median(dx)),
            dy_median=float(np.median(dy)),
            offset_median=offset,
            offset_mad=float(np.median(np.abs(resid - offset))),
            n_frame_stars=int(len(fx)),
        )

    if not groups:
        raise GaugeError(
            f"no (visit, detector) group reached {min_matched} matches against "
            f"{i2d_path}; check --filter/--detector/--pupil and --exposure "
            f"{reference_exposure}"
        )

    result = GaugeResult(
        groups=groups,
        i2d_scale_arcsec=scale,
        i2d_shape=(ny, nx),
        match_arcsec=match_arcsec,
        n_i2d_stars=int(len(ix)),
        skipped=skipped,
    )
    result.visits = summarise_visits(groups)
    result.anchor_visit = min(
        result.visits, key=lambda v: result.visits[v].offset_median
    )
    return result


def summarise_visits(
    groups: dict[str, GaugeGroup],
) -> dict[int, GaugeVisit]:
    """Roll per-group offsets up to one row per visit (median of medians).

    The median is taken across detectors rather than pooling stars, so a
    detector with more matches cannot dominate the visit.  ``spread`` is the
    max-minus-min of the group's ``offset_median`` and ``internal_rms`` is the
    rms of the per-detector vectors about the visit centre: the first is how
    much the detectors disagree, the second how tightly they cluster.  Both are
    diagnostics for whether the visit median can be trusted.
    """
    visits: dict[int, GaugeVisit] = {}
    for visit in sorted({g.visit for g in groups.values()}):
        members = [g for g in groups.values() if g.visit == visit]
        dx = float(np.median([g.dx_median for g in members]))
        dy = float(np.median([g.dy_median for g in members]))
        offsets = [g.offset_median for g in members]
        dev = np.hypot(
            np.array([g.dx_median for g in members]) - dx,
            np.array([g.dy_median for g in members]) - dy,
        )
        visits[visit] = GaugeVisit(
            visit=visit,
            dx_median=dx,
            dy_median=dy,
            offset_median=float(np.hypot(dx, dy)),
            spread=float(max(offsets) - min(offsets)),
            n_det=len(members),
            internal_rms=float(np.sqrt(np.mean(dev**2))),
        )
    return visits


def suspect_groups(result: GaugeResult) -> list[str]:
    """Groups measured too close to the edge of the search window to trust.

    ``n_matched`` and scatter both fail to detect an offset wider than the
    radius, and no cheap statistic does -- see
    :data:`SUSPECT_RADIUS_FRACTION`.  What can be checked is whether the
    *median* offset has grown to a sizeable fraction of the radius, which is
    the regime where a real translation and a lucky neighbour pairing are
    indistinguishable.  All shipped filters sit near 10%.
    """
    if result.i2d_scale_arcsec <= 0:
        return []
    radius_px = result.match_arcsec / result.i2d_scale_arcsec
    return sorted(
        key
        for key, g in result.groups.items()
        if g.offset_median > SUSPECT_RADIUS_FRACTION * radius_px
    )


def format_gauge_report(result: GaugeResult, *, filter_label: str = "") -> str:
    """Human-readable report, in the layout the recorded step-2 logs use."""
    lines: list[str] = []
    label = f"filter={filter_label}" if filter_label else ""
    lines.append(f"step 2 gauge probe {label}".rstrip())
    ny, nx = result.i2d_shape
    lines.append(
        f"i2d shape {ny} x {nx}, measured scale {result.i2d_scale_arcsec:.6f} "
        f"arcsec/px"
    )
    lines.append(
        f"match radius {result.match_arcsec} arcsec = "
        f"{result.match_arcsec / result.i2d_scale_arcsec:.2f} i2d px"
    )
    lines.append(f"  {result.n_i2d_stars} stars in the i2d")
    for s in result.skipped:
        lines.append(f"  SKIPPED {s}")
    lines.append("")
    lines.append("per (visit x detector) raw frame WCS vs official i2d, in i2d px:")
    lines.append("  visit detector   n_match  rate    dx_med    dy_med   |offset|    scatter")
    suspect = suspect_groups(result)
    for key in sorted(result.groups, key=lambda k: (result.groups[k].visit,
                                                    result.groups[k].detector)):
        g = result.groups[key]
        lines.append(
            f"  v{g.visit}   {g.detector:<10} {g.n_matched:4d}  "
            f"{g.match_rate:5.0%}  {g.dx_median:+8.4f}  {g.dy_median:+8.4f}  "
            f"{g.offset_median:7.4f}  {g.offset_mad:8.4f}"
        )
    if suspect:
        radius_px = result.match_arcsec / result.i2d_scale_arcsec
        lines.append("")
        lines.append(
            f"WARNING: {len(suspect)} group(s) sit within "
            f"{SUSPECT_RADIUS_FRACTION:.0%} of the {radius_px:.2f} px match "
            f"radius, so a real offset and a lucky neighbour pairing are not "
            f"distinguishable here:"
        )
        for key in suspect:
            g = result.groups[key]
            lines.append(
                f"  {key}: |offset| {g.offset_median:.4f} px of a "
                f"{radius_px:.2f} px radius, {g.n_matched} matched. Re-check "
                f"--match-radius-arcsec and the match rate before using this."
            )
    lines.append("")
    n_det = max((v.n_det for v in result.visits.values()), default=0)
    lines.append(
        f"per-visit summary across the {n_det} detector(s) (median of medians):"
    )
    lines.append("  visit   dx_med   dy_med  |offset|   spread(max-min)  n_det")
    for visit in sorted(result.visits):
        v = result.visits[visit]
        lines.append(
            f"  v{visit}   {v.dx_median:+8.4f}  {v.dy_median:+8.4f}  "
            f"{v.offset_median:7.4f}      {v.spread:6.4f}        "
            f"{v.n_det}   (internal {v.internal_rms:.4f})"
        )
    if result.anchor_visit is not None:
        anchor = result.visits[result.anchor_visit]
        lines.append("")
        lines.append(
            f"anchor by per-visit median: visit {result.anchor_visit} "
            f"({anchor.offset_median:.4f} i2d px "
            f"= {px_to_arcsec(anchor.offset_median, result.i2d_scale_arcsec, I2D_PX):.5f} arcsec)"
        )
        lines.append("")
        lines.append(
            "visit-to-anchor separation, i.e. the raw inter-visit pointing error "
            "this mosaic must absorb:"
        )
        for visit in sorted(result.visits):
            sep = result.separation_from_anchor(visit)
            mark = "   <- anchor" if visit == result.anchor_visit else ""
            lines.append(
                f"  v{visit}: {sep:.4f} i2d px "
                f"= {px_to_arcsec(sep, result.i2d_scale_arcsec, I2D_PX):.5f} arcsec{mark}"
            )
        if len(result.visits) > 1:
            worst = max(
                (v for v in result.visits if v != result.anchor_visit),
                key=result.separation_from_anchor,
                default=None,
            )
            if worst is not None:
                lines.append(
                    f"  worst visit: v{worst} at "
                    f"{result.separation_from_anchor(worst):.4f} px raw"
                )
    return "\n".join(lines)


def write_gauge_outputs(
    result: GaugeResult,
    *,
    json_path: str | Path,
    log_path: str | Path | None = None,
    filter_label: str = "",
) -> None:
    """Write the flat per-group JSON and, optionally, the text report.

    The JSON is the flat ``{"v1_nrca1": {...}}`` form used by the existing
    ``*_step2_detectors.json`` artifacts, so it is a peer of them rather than a
    competing schema.
    """
    json_path = Path(json_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(result.as_groups_dict(), indent=2) + "\n", encoding="utf-8"
    )
    if log_path is not None:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            format_gauge_report(result, filter_label=filter_label) + "\n",
            encoding="utf-8",
        )


def filter_label(exposures: Sequence[CalExposure]) -> str:
    """``F444W;CLEAR``-style label for the exposures actually selected."""
    filters = sorted({e.filter_name.strip().upper() for e in exposures})
    if not filters:
        return ""
    if len(filters) > 1:
        return ";".join(filters)
    pupils = sorted({e.pupil.strip().upper() for e in exposures})
    if len(pupils) == 1 and pupils[0] not in ("", "CLEAR"):
        return f"{filters[0]};{pupils[0]}"
    return filters[0]
