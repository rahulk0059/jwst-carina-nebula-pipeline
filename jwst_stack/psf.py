"""Per-channel stellar FWHM, measured on the finished mosaics.

This is the Phase 1 measurement, written as re-runnable code rather than
transcribed from a throwaway script.  It answers one question per channel:
**how wide is the PSF on the 0.031 arcsec colour grid?**  That number is the
honest limit on what the composite can show, it is what says whether the three
RGB channels can be displayed at equal weight (a 2.2x FWHM spread between
F090W and F470N does not mean F470N is 2.2x sharper, it means its stars are
2.2x more spread out on the grid), and it belongs in ``color_channels.json``
next to the backgrounds so a reader can see the two calibration facts together.

Method, and why each piece is there:

* **Tiled sampling over the field, not a crop.**  :func:`jwst_stack.starcat.
  detect_tiled` runs the same DAOStarFinder the rest of the project uses over a
  spread of well-covered tiles, so the sample is drawn from across the mosaic
  rather than from one detector's corner.  Detection runs on a halo-expanded
  tile and keeps only interior centroids, which makes the catalogue independent
  of the tiling.
* **A 2-D rotated Gaussian per star, plus a local pedestal.**  The pedestal is
  fitted rather than subtracted, so the answer is self-referenced exactly as
  the background policy in :mod:`jwst_stack.color` is: no external level, no
  i2d, nothing that could differ between the pass that was measured and a
  rerun.  The patch half-width grows as 4 sigma of the expected width, so a
  5.5 px FWHM is not fitted inside a box too small to contain it.
* **A second detection pass at the measured width.**  DAOStarFinder's ``fwhm``
  argument sets the kernel it convolves with, so using a fixed 3 px guess
  under-detects the 5.5 px channels.  Pass 1 measures, pass 2 re-detects with
  that measurement, and the difference between the passes is itself reported -
  agreement is evidence the detection is not setting the answer.
* **A median over separated stars.**  Stars are accepted brightest-first and
  rejected within a separation of each other, so blends and double cores are
  dropped rather than averaged in, and the aggregate is a median so a few
  outliers cannot move it.  The reported FWHM is the geometric mean of the two
  fitted axes, which is the natural scalar for an elliptical PSF; the
  axis ratio is recorded separately as a quality diagnostic.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit

from jwst_stack.io import update_json_section
from jwst_stack.starcat import (
    DEFAULT_DETECT_FWHM_PX,
    DEFAULT_THRESHOLD_SIGMA,
    CoverageProbe,
    brightest_separated,
    covered_tiles,
    detect_tiled,
)
from jwst_stack.units import GRID_PX, px_to_arcsec

#: Full width of a Gaussian at half maximum, in units of its sigma.
FWHM_PER_SIGMA = 2.0 * math.sqrt(2.0 * math.log(2.0))

#: Smallest patch half-width, in pixels, for a narrow PSF.
MIN_HALF_WIDTH_PX = 6

#: Patch half-width in units of sigma, so a wide PSF is fitted on a wide patch.
PATCH_SIGMA_FACTOR = 4.0

#: Two stars closer than this are treated as one blended source and only the
#: brighter one is fitted.
DEFAULT_SEPARATION_PX = 12.0

#: Stars to fit per channel.  Far more than the number that survive sigma
#: clipping, so the median is not decided by a handful of stars.
DEFAULT_TARGET_STARS = 200

#: A channel with fewer fitted stars than this is not reported as measured.
DEFAULT_MIN_STARS = 8

#: Clip for the per-star aggregate, in sigma of the sample's own spread.
DEFAULT_STAR_CLIP_SIGMA = 3.0

#: Sigma outside which a single star's fitted width is rejected outright.
MAX_REASONABLE_FWHM_PX = 30.0

#: Halo floor for the FWHM pass.  The patch size normally sets the real halo,
#: since a star needs the whole patch inside the expanded tile.
DEFAULT_HALO_PX_FOR_PATCH = 24


def patch_half_width(fwhm_px: float) -> int:
    """Half-width in pixels of a patch that comfortably contains *fwhm_px*."""
    sigma = max(float(fwhm_px), 1e-6) / FWHM_PER_SIGMA
    return int(max(MIN_HALF_WIDTH_PX, math.ceil(PATCH_SIGMA_FACTOR * sigma)))


def gaussian_2d(
    xy: np.ndarray,
    pedestal: float,
    amp: float,
    x0: float,
    y0: float,
    sig_x: float,
    sig_y: float,
    theta: float,
) -> np.ndarray:
    """Rotated 2-D Gaussian on a flat pedestal, evaluated at stacked ``(x, y)``.

    *xy* is the ``np.stack([xx, yy])`` pair, both raveled - ``curve_fit`` wants
    one leading coordinate axis, and passing a 2-D meshgrid makes it treat each
    pixel as an independent parameter.
    """
    ct = math.cos(theta)
    st = math.sin(theta)
    dx = xy[0] - x0
    dy = xy[1] - y0
    u = dx * ct + dy * st
    v = -dx * st + dy * ct
    return pedestal + amp * np.exp(-0.5 * ((u / sig_x) ** 2 + (v / sig_y) ** 2))


@dataclass
class StarFit:
    """One fitted star."""

    x: float
    y: float
    flux: float
    fwhm_px: float
    fwhm_major_px: float
    fwhm_minor_px: float
    axis_ratio: float
    amp: float
    pedestal: float
    rms_residual: float


def fit_star_patch(
    patch: np.ndarray,
    cx: float,
    cy: float,
    *,
    fwhm_hint_px: float = DEFAULT_DETECT_FWHM_PX,
    x: float = 0.0,
    y: float = 0.0,
    flux: float = float("nan"),
) -> StarFit | None:
    """Fit one star in *patch*, centred at pixel (*cx*, *cy*).

    *patch* must be entirely finite.  Returns ``None`` when the fit does not
    converge or lands on an unphysical width, so a bad star is dropped rather
    than silently entering the median.
    """
    arr = np.asarray(patch, dtype=float)
    if arr.ndim != 2 or not np.all(np.isfinite(arr)):
        return None
    ny, nx = arr.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    xy = np.stack([xx.ravel(), yy.ravel()], axis=0)
    flat = arr.ravel()

    sigma_hint = max(float(fwhm_hint_px), 1e-6) / FWHM_PER_SIGMA
    pedestal0 = float(np.median(flat))
    amp0 = float(flat.max() - pedestal0)
    if not np.isfinite(amp0) or amp0 <= 0.0:
        return None
    p0 = [pedestal0, amp0, float(cx), float(cy), sigma_hint, sigma_hint, 0.0]
    span = float(min(nx, ny)) / 2.0
    lower = [-np.inf, 0.0, cx - span, cy - span, 0.2, 0.2, -math.pi / 2]
    upper = [np.inf, np.inf, cx + span, cy + span, 50.0, 50.0, math.pi / 2]
    try:
        popt, _ = curve_fit(gaussian_2d, xy, flat, p0=p0, bounds=(lower, upper), maxfev=4000)
    except (RuntimeError, ValueError, TypeError):
        return None
    pedestal, amp, x0, y0, sig_x, sig_y, theta = (float(v) for v in popt)
    if not all(np.isfinite(v) for v in (pedestal, amp, x0, y0, sig_x, sig_y, theta)):
        return None
    if amp <= 0.0 or sig_x <= 0.0 or sig_y <= 0.0:
        return None
    major = FWHM_PER_SIGMA * max(sig_x, sig_y)
    minor = FWHM_PER_SIGMA * min(sig_x, sig_y)
    if not (0.0 < minor <= MAX_REASONABLE_FWHM_PX and major <= MAX_REASONABLE_FWHM_PX):
        return None
    residual = gaussian_2d(xy, pedestal, amp, x0, y0, sig_x, sig_y, theta) - flat
    rms = float(np.sqrt(np.mean(residual**2)))
    return StarFit(
        x=x,
        y=y,
        flux=flux,
        fwhm_px=math.sqrt(major * minor),
        fwhm_major_px=major,
        fwhm_minor_px=minor,
        axis_ratio=major / minor if minor > 0 else float("nan"),
        amp=amp,
        pedestal=pedestal,
        rms_residual=rms,
    )


@dataclass
class FwhmResult:
    """One channel's measured PSF width."""

    channel: str
    mosaic: str
    n_tiles: int
    n_candidates: int
    n_fitted: int
    n_used: int
    detect_fwhm_px: float
    fwhm_px: float = float("nan")
    fwhm_arcsec: float = float("nan")
    fwhm_mad_px: float = float("nan")
    axis_ratio_median: float = float("nan")
    grid_scale_arcsec: float = float("nan")
    per_star_px: list[float] = field(default_factory=list)

    @property
    def measured(self) -> bool:
        return self.n_used >= DEFAULT_MIN_STARS

    def as_dict(self) -> dict:
        return {
            "channel": self.channel,
            "mosaic": self.mosaic,
            "grid_scale_arcsec_per_px": self.grid_scale_arcsec,
            "pixel_basis": GRID_PX,
            "detect_fwhm_px": self.detect_fwhm_px,
            "fwhm_px": self.fwhm_px,
            "fwhm_arcsec": self.fwhm_arcsec,
            "fwhm_mad_px": self.fwhm_mad_px,
            "axis_ratio_median": self.axis_ratio_median,
            "n_tiles": self.n_tiles,
            "n_candidates": self.n_candidates,
            "n_fitted": self.n_fitted,
            "n_used": self.n_used,
            "measured": self.measured,
            "per_star_fwhm_px": [float(v) for v in self.per_star_px],
        }


def _fit_stars_on_tiles(
    image: np.ndarray,
    shape: tuple[int, int],
    tiles: list[tuple[int, int, int, int]],
    *,
    fwhm_hint_px: float,
    target_stars: int,
    separation_px: float,
    threshold_sigma: float,
    half_width: int,
    per_tile_stars: int,
    max_detected: int,
    progress=None,
) -> list[StarFit]:
    """Cut, fit and clip patches for every separated star in *tiles*."""
    x, y, flux = detect_tiled(
        image,
        tiles,
        halo_px=max(DEFAULT_HALO_PX_FOR_PATCH, half_width + 4),
        keep_margin=float(half_width) + 1.0,
        fwhm_px=fwhm_hint_px,
        threshold_sigma=threshold_sigma,
        per_tile_stars=per_tile_stars,
        max_stars=max_detected,
        progress=progress,
    )
    picks = brightest_separated(x, y, flux, target_stars, separation_px)
    fits: list[StarFit] = []
    arr = np.asarray(image)
    for idx in picks:
        px = int(round(float(x[idx])))
        py = int(round(float(y[idx])))
        y0, y1 = py - half_width, py + half_width + 1
        x0, x1 = px - half_width, px + half_width + 1
        patch = np.asarray(arr[y0:y1, x0:x1], dtype=float)
        if patch.shape != (2 * half_width + 1, 2 * half_width + 1):
            continue
        fit = fit_star_patch(
            patch,
            float(x[idx]) - x0,
            float(y[idx]) - y0,
            fwhm_hint_px=fwhm_hint_px,
            x=float(x[idx]),
            y=float(y[idx]),
            flux=float(flux[idx]),
        )
        if fit is not None:
            fits.append(fit)
    return fits


def measure_fwhm(
    image: np.ndarray,
    coverage: CoverageProbe | None,
    shape: tuple[int, int],
    *,
    channel: str = "",
    mosaic: str = "",
    grid_scale_arcsec: float = float("nan"),
    n_tiles: int | None = 24,
    tile_px: int = 512,
    separation_px: float = DEFAULT_SEPARATION_PX,
    target_stars: int = DEFAULT_TARGET_STARS,
    threshold_sigma: float = DEFAULT_THRESHOLD_SIGMA,
    clip_sigma: float = DEFAULT_STAR_CLIP_SIGMA,
    min_coverage: float = 0.9,
    fwhm_hint_px: float = DEFAULT_DETECT_FWHM_PX,
    per_tile_stars: int = 200,
    max_detected: int = 600,
    progress=None,
) -> FwhmResult:
    """Measure the stellar FWHM of one channel on the mosaic grid.

    *image* is a memmapped view and *coverage* is a probe from
    :func:`jwst_stack.starcat.coverage_probe`; only the sampled tiles are read.
    Pass *n_tiles* as ``None`` for a full-footprint sample (slower, but it
    removes the sampling question entirely).
    """
    tiles = covered_tiles(shape, coverage, tile_px, n_tiles=n_tiles, min_coverage=min_coverage)
    half_width = patch_half_width(fwhm_hint_px)
    fits = _fit_stars_on_tiles(
        image,
        shape,
        tiles,
        fwhm_hint_px=fwhm_hint_px,
        target_stars=target_stars,
        separation_px=separation_px,
        threshold_sigma=threshold_sigma,
        half_width=half_width,
        per_tile_stars=per_tile_stars,
        max_detected=max_detected,
        progress=progress,
    )
    result = FwhmResult(
        channel=channel,
        mosaic=mosaic,
        n_tiles=len(tiles),
        n_candidates=len(fits),
        n_fitted=len(fits),
        n_used=0,
        detect_fwhm_px=float(fwhm_hint_px),
        grid_scale_arcsec=float(grid_scale_arcsec),
    )
    if not fits:
        return result
    values = np.asarray([f.fwhm_px for f in fits], dtype=float)
    if values.size >= 4:
        clipped = sigma_clip(values, sigma=clip_sigma, maxiters=3, masked=True)
        keep = ~np.ma.getmaskarray(clipped)
    else:
        keep = np.ones(values.shape, dtype=bool)
    used = values[keep]
    if used.size == 0:
        return result
    median = float(np.median(used))
    result.fwhm_px = median
    result.fwhm_mad_px = 1.4826 * float(np.median(np.abs(used - median)))
    result.n_used = int(used.size)
    result.axis_ratio_median = float(
        np.median([f.axis_ratio for f, k in zip(fits, keep) if k])
    )
    if np.isfinite(grid_scale_arcsec):
        result.fwhm_arcsec = px_to_arcsec(median, grid_scale_arcsec, GRID_PX)
    result.per_star_px = [float(v) for v in used]
    return result


def measure_channels_fwhm(
    mosaics: dict[str, str],
    *,
    grid_scale_arcsec: float = float("nan"),
    refine: bool = True,
    progress=None,
    **kwargs,
) -> tuple[dict[str, FwhmResult], dict[str, FwhmResult]]:
    """Measure every channel, then re-detecting at the measured width.

    With *refine*, each channel is measured twice: once with the fixed 3 px
    detection hint, then again with the first pass's answer as the hint.  The
    second pass is returned, because DAOStarFinder's kernel width otherwise
    biases the wide channels low.

    Returns ``(results, first_pass)``.  Both are returned rather than only the
    final one so the refinement can be reported: the gap between the two passes
    is the size of the detection bias, which is evidence about the measurement
    rather than something to hide.
    """
    from jwst_stack.starcat import open_mosaic

    base_hint = float(kwargs.pop("fwhm_hint_px", DEFAULT_DETECT_FWHM_PX))
    first: dict[str, FwhmResult] = {}
    for name, path in mosaics.items():
        data = open_mosaic(path)
        first[name] = measure_fwhm(
            data.sci,
            data.coverage,
            data.shape,
            channel=name,
            mosaic=str(path),
            grid_scale_arcsec=grid_scale_arcsec,
            fwhm_hint_px=base_hint,
            progress=progress,
            **kwargs,
        )
    if not refine:
        return first, first
    refined: dict[str, FwhmResult] = {}
    for name, path in mosaics.items():
        hint = first[name].fwhm_px
        if not np.isfinite(hint) or hint <= 0:
            refined[name] = first[name]
            continue
        data = open_mosaic(path)
        refined[name] = measure_fwhm(
            data.sci,
            data.coverage,
            data.shape,
            channel=name,
            mosaic=str(path),
            grid_scale_arcsec=grid_scale_arcsec,
            fwhm_hint_px=float(hint),
            progress=progress,
            **kwargs,
        )
    return refined, first


def write_fwhm_json(
    path: str | Path,
    results: dict[str, FwhmResult],
    *,
    first_pass: dict[str, FwhmResult] | None = None,
    grid_path: str | Path | None = None,
    grid_shape: tuple[int, int] | None = None,
    extra: dict | None = None,
) -> Path:
    """Write the ``fwhm`` section of ``color_channels.json``.

    Both passes are recorded when available, so the refinement is auditable
    rather than implied.  Existing top-level sections are preserved: this file
    is also where the backgrounds and the stretch parameters end up.
    """
    payload = {
        "method": (
            "2-D rotated Gaussian + local pedestal per star, on separated stars "
            "sampled across the field; median over sigma-clipped stars; "
            "DAOStarFinder re-run at the pass-1 width"
        ),
        "grid": {
            "path": str(grid_path) if grid_path else None,
            "scale_arcsec_per_px": float(grid_scale_arcsec_or_nan(results)),
            "shape": list(grid_shape) if grid_shape else None,
        },
        "channels": {name: result.as_dict() for name, result in results.items()},
    }
    if first_pass is not None:
        payload["pass1_detect_hint_3px"] = {
            name: {
                "fwhm_px": r.fwhm_px,
                "fwhm_arcsec": r.fwhm_arcsec,
                "n_used": r.n_used,
            }
            for name, r in first_pass.items()
        }
    if extra:
        payload.update(extra)
    return update_json_section(path, "fwhm", payload)


def grid_scale_arcsec_or_nan(results: dict[str, FwhmResult]) -> float:
    """The grid scale, taken from whichever channel recorded it."""
    for result in results.values():
        if np.isfinite(result.grid_scale_arcsec):
            return float(result.grid_scale_arcsec)
    return float("nan")


def format_fwhm_summary(results: dict[str, FwhmResult]) -> str:
    """Terminal table, basis-tagged like every other pixel number here."""
    lines = ["stellar FWHM on the colour grid (grid_px)"]
    header = (
        f"  {'channel':<14} {'FWHM px':>9} {'FWHM arcsec':>12} {'MAD px':>7} "
        f"{'axis':>6} {'stars':>7} {'tiles':>6}"
    )
    lines.append(header)
    for name, r in results.items():
        if not r.measured:
            lines.append(f"  {name:<14} {'not measured':>9} ({r.n_fitted} stars fitted)")
            continue
        lines.append(
            f"  {name:<14} {r.fwhm_px:9.2f} {r.fwhm_arcsec:12.5f} "
            f"{r.fwhm_mad_px:7.3f} {r.axis_ratio_median:6.3f} {r.n_used:7d} {r.n_tiles:6d}"
        )
    return "\n".join(lines)
