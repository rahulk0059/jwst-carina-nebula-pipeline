"""Per-channel background subtraction and a colour-preserving display stretch.

The six mosaics are co-located on the shared 0.031 arcsec grid and can be shown
as one RGB (F090W / F200W / F444W;CLEAR) plus three separate feature panels
(F187N, F335M, F444W;F470N).  This module holds the two decisions that determine
whether that composite is honest:

1. **Background is self-referenced, per channel.**  Each channel has its own
   additive pedestal estimated from its own pixels and subtracted from its own
   pixels.  No channel is scaled or offset to match another, and the official
   i2d is never consulted.  That is what makes the mosaic-vs-i2d background
   differences - F187N at +4.8% and F335M at +4.5% of the official sky, the
   largest in the set - irrelevant to colour balance.  Those percentages
   describe *our document minus the official one*; subtracting our own measured
   level removes our pedestal whatever its size, and the official pedestal never
   enters.  Concretely, :func:`source_free_level` is translation-covariant
   (median and MAD both shift with the data), so for any constant ``c``

       ``level(x + c) == level(x) + c``   and therefore   ``(x + c) - level(x + c) == x - level(x)``

   exactly.  An additive pedestal in any channel, of any size, cannot move the
   displayed colour.  ``test_additive_pedestal_does_not_change_colour`` pins it.

2. **The stretch is shared, not per channel.**  After subtraction the clipping
   limits are computed once, from all channels pooled, and the same ``lo``/``hi``
   and the same softening ``a`` are applied to every channel.  Per-channel
   min-max normalisation - which is what :func:`plotting.asinh_stretch` does -
   would destroy exactly the relative brightness that carries the colour, by
   silently rescaling a faint channel and a bright one to the same [0, 1] range.

Background here is the *typical scene level*, estimated robustly so that stars
and cosmic rays do not drag it up.  It is not a dark-sky-only estimate: in a
nebula-filled field the bulk level *is* the nebula, and centring the asinh on it
is what puts the mid-tones in the middle of the display range.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np

from jwst_stack.mosaic import iter_tiles
from jwst_stack.units import GRID_PX, px_label, px_to_arcsec

#: Tile edge for the robust level estimate.  Only local statistics are taken per
#: tile, so this bounds memory, not the answer.
DEFAULT_TILE_PX = 512

#: Clip distance in sigma (1.4826 x MAD, the convention used in validation.py).
DEFAULT_CLIP_SIGMA = 3.0

#: A tile with fewer finite pixels than this contributes nothing.
DEFAULT_MIN_PIXELS = 16

#: Shared percentile clip, applied to every channel with the same limits.
DEFAULT_LO_PCT = 1.0
DEFAULT_HI_PCT = 99.5

#: Softening ``a`` as a fixed fraction of the shared display range.
DEFAULT_SOFTENING_FRACTION = 0.1

#: Cap on the pixels pooled per channel when finding the shared limits.
DEFAULT_SAMPLE_CAP = 2_000_000


def _sigma_mad(values: np.ndarray) -> float:
    """1.4826 x MAD, matching :func:`jwst_stack.validation._sigma_mad`."""
    med = float(np.median(values))
    return 1.4826 * float(np.median(np.abs(values - med)))


def source_free_level(
    image: np.ndarray,
    *,
    tile_px: int = DEFAULT_TILE_PX,
    clip_sigma: float = DEFAULT_CLIP_SIGMA,
    min_pixels: int = DEFAULT_MIN_PIXELS,
) -> float:
    """Robust typical-scene level of *image*, ignoring sources and holes.

    Each tile is sigma-clipped about its own median (so stars, cosmic rays and
    the bright tail do not bias it) and contributes one clipped median; the
    result is the median of those per-tile medians.  Per-tile then median, rather
    than one pooling of every kept pixel, so a large bright region cannot outvote
    the rest of the field.

    Reads tile by tile, so a memmapped mosaic is streamed rather than loaded.
    Returns NaN when the image has no usable pixels.
    """
    arr = np.asarray(image)
    levels: list[float] = []
    for y0, y1, x0, x1 in iter_tiles(arr.shape, tile_px):
        tile = np.asarray(arr[y0:y1, x0:x1], dtype=float)
        values = tile[np.isfinite(tile)]
        if values.size < min_pixels:
            continue
        median = float(np.median(values))
        sigma = _sigma_mad(values)
        if sigma > 0.0:
            kept = values[np.abs(values - median) <= clip_sigma * sigma]
            if kept.size >= min_pixels:
                values = kept
                median = float(np.median(values))
        levels.append(median)
    if not levels:
        return float("nan")
    return float(np.median(levels))


def subtract_background(image: np.ndarray, background: float) -> np.ndarray:
    """``image - background``, preserving dtype and NaN holes."""
    return np.asarray(image) - background


def shared_limits(
    images: Mapping[str, np.ndarray] | list[np.ndarray],
    *,
    lo_pct: float = DEFAULT_LO_PCT,
    hi_pct: float = DEFAULT_HI_PCT,
    sample_cap: int = DEFAULT_SAMPLE_CAP,
) -> tuple[float, float]:
    """Percentile limits pooled across *images*, so every channel shares them.

    The inputs are expected to be background-subtracted already.  Each channel
    contributes at most ``sample_cap`` evenly spaced finite samples, which keeps
    the helper usable on a memmapped mosaic without materialising it.
    """
    if isinstance(images, Mapping):
        arrays: list[np.ndarray] = list(images.values())
    else:
        arrays = list(images)
    samples: list[np.ndarray] = []
    for image in arrays:
        flat = np.asarray(image).ravel()
        if flat.size == 0:
            continue
        step = max(1, flat.size // sample_cap)
        sub = np.asarray(flat[::step], dtype=float)
        sub = sub[np.isfinite(sub)]
        if sub.size:
            samples.append(sub)
    if not samples:
        return (float("nan"), float("nan"))
    lo, hi = np.percentile(np.concatenate(samples), [lo_pct, hi_pct])
    return (float(lo), float(hi))


def default_softening(
    lo: float, hi: float, fraction: float = DEFAULT_SOFTENING_FRACTION
) -> float:
    """Softening ``a`` as a fixed fraction of the shared display range.

    Fixed because it is shared: the same ``a`` on every channel is what keeps the
    stretch linear in relative brightness over the faint end.
    """
    return float(fraction) * (hi - lo)


def background_centred_asinh(
    image: np.ndarray,
    *,
    background: float,
    a: float,
    lo: float,
    hi: float,
) -> np.ndarray:
    """Asinh stretch centred on the measured *background*, clipped by *lo*/*hi*.

    ``background`` maps to ``arcsinh(0) = 0`` and the shared ``lo``/``hi`` map to
    0 and 1, so the clip range means the same thing in every channel.  NaN holes
    stay NaN; finite values are clipped to [0, 1].
    """
    if not np.isfinite(a) or a <= 0.0:
        raise ValueError(f"softening must be positive and finite, got {a!r}")
    if not np.isfinite(hi) or not np.isfinite(lo) or hi <= lo:
        raise ValueError(f"display limits must satisfy lo < hi, got {lo!r}, {hi!r}")
    arr = np.asarray(image, dtype=float)
    y = np.arcsinh((arr - background) / a)
    y_lo = np.arcsinh((lo - background) / a)
    y_hi = np.arcsinh((hi - background) / a)
    out = (y - y_lo) / (y_hi - y_lo)
    return np.clip(out, 0.0, 1.0)


@dataclass
class StretchResult:
    """Background-subtracted, stretched channels plus the parameters used."""

    arrays: dict[str, np.ndarray]
    backgrounds: dict[str, float]
    lo: float
    hi: float
    a: float
    softening_fraction: float

    def as_dict(self) -> dict:
        return {
            "backgrounds": {k: float(v) for k, v in self.backgrounds.items()},
            "lo": self.lo,
            "hi": self.hi,
            "softening": self.a,
            "softening_fraction": self.softening_fraction,
        }


def stretch_channels(
    images: Mapping[str, np.ndarray],
    *,
    lo_pct: float = DEFAULT_LO_PCT,
    hi_pct: float = DEFAULT_HI_PCT,
    softening_fraction: float = DEFAULT_SOFTENING_FRACTION,
    tile_px: int = DEFAULT_TILE_PX,
    clip_sigma: float = DEFAULT_CLIP_SIGMA,
) -> StretchResult:
    """Apply the shared-stretch policy to a set of named channels.

    Subtracts each channel's own :func:`source_free_level`, then stretches all of
    them with one shared pair of limits and one shared softening.  In-memory: it
    holds the subtracted arrays, so Phase 4 streams this over tiles for the
    0.031 grid; here it is the readable statement of the policy and the target of
    the invariance tests.
    """
    backgrounds = {
        name: source_free_level(image, tile_px=tile_px, clip_sigma=clip_sigma)
        for name, image in images.items()
    }
    subtracted = {
        name: subtract_background(image, backgrounds[name])
        for name, image in images.items()
    }
    lo, hi = shared_limits(subtracted, lo_pct=lo_pct, hi_pct=hi_pct)
    a = default_softening(lo, hi, softening_fraction)
    arrays = {
        name: background_centred_asinh(image, background=0.0, a=a, lo=lo, hi=hi)
        for name, image in subtracted.items()
    }
    return StretchResult(
        arrays=arrays,
        backgrounds=backgrounds,
        lo=lo,
        hi=hi,
        a=a,
        softening_fraction=softening_fraction,
    )


# --------------------------------------------------------------------------
# the six channels
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ChannelSpec:
    """One display channel and the mosaic it is built from.

    ``role`` is ``R``/``G``/``B`` for the three that go into the composite and
    ``feature`` for the narrowbands that get their own panel.  Keeping the
    feature channels in the same table is deliberate: they are part of the same
    six-product measurement, and leaving them out of the RGB is a property of
    this table rather than a rule the builder has to remember.
    """

    name: str
    filter_name: str
    pupil: str
    role: str
    mosaic: str

    @property
    def bandpass(self) -> str:
        """``FILTER`` and ``PUPIL`` as a product name, e.g. ``F444W;F470N``."""
        if self.pupil and self.pupil.upper() != "CLEAR":
            return f"{self.filter_name};{self.pupil}"
        return self.filter_name

    @property
    def is_feature(self) -> bool:
        return self.role == "feature"

    def mosaic_path(self, outdir: str | Path) -> Path:
        return Path(outdir) / self.mosaic

    def validation_json(self, outdir: str | Path) -> Path:
        """The ``compare`` JSON for this channel, if it has been run.

        F335M, F444W and F444W;F470N have both a native-grid and a 0.031-grid
        mosaic, so the file is found by globbing rather than by constructing
        the name.  The colour build consumes the 0.031-grid product, so the
        grid marker is preferred where both exist.
        """
        stem = self.mosaic[: -len(".fits")]
        candidates = sorted(Path(outdir).glob(f"{stem}_vs_*.json"))
        if not candidates:
            return Path(outdir) / f"{stem}_vs_<missing>.json"
        grid_markers = [p for p in candidates if "0031grid" in p.name]
        return (grid_markers or candidates)[0]


#: The six channels, in display order: RGB first, then the feature panels.
CHANNELS: tuple[ChannelSpec, ...] = (
    ChannelSpec("F444W", "F444W", "CLEAR", "R", "f444w_all_detectors_mosaic_0031grid.fits"),
    ChannelSpec("F200W", "F200W", "CLEAR", "G", "f200w_all_detectors_mosaic.fits"),
    ChannelSpec("F090W", "F090W", "CLEAR", "B", "f090w_all_detectors_mosaic.fits"),
    ChannelSpec("F187N", "F187N", "CLEAR", "feature", "f187n_all_detectors_mosaic.fits"),
    ChannelSpec("F335M", "F335M", "CLEAR", "feature", "f335m_all_detectors_mosaic_0031grid.fits"),
    ChannelSpec(
        "F444W;F470N",
        "F444W",
        "F470N",
        "feature",
        "f444w_f470n_all_detectors_mosaic_0031grid.fits",
    ),
)

#: The registration reference for the colour build.
COLOR_REFERENCE = "F200W"

_CHANNEL_BY_NAME: dict[str, ChannelSpec] = {spec.name: spec for spec in CHANNELS}

#: Match radius for the cross-filter star match, in colour-grid pixels.  4 px at
#: 0.031 arcsec/px is 0.124 arcsec, comfortably inside the ~2.5 px FWHM of the
#: narrowest channel, so genuine counterparts match and unrelated neighbours do
#: not.
DEFAULT_MATCH_RADIUS_PX = 4.0

#: A median cross-filter shift below this is recorded but not worth correcting.
#: It is a quarter pixel, or 0.0078 arcsec on this grid - below the ~0.05 px
#: comparison floor the project has measured against these i2d products, so
#: shifting by it would be fitting noise.
DEFAULT_NEGLIGIBLE_SHIFT_PX = 0.25

#: Detection settings for the cross-filter match.  A full-footprint catalogue
#: is used: a sampled one would be a poorer match than the Phase 2 run it
#: replaces, and the star count is the thing that makes the median trustworthy.
DEFAULT_REGISTER_TILE_PX = 1024
DEFAULT_REGISTER_HALO_PX = 24
DEFAULT_REGISTER_DETECT_FWHM_PX = 3.0


def channel_table(outdir: str | Path) -> dict[str, ChannelSpec]:
    """The six channels keyed by display name."""
    return {spec.name: spec for spec in CHANNELS}


def rgb_channels() -> list[ChannelSpec]:
    """The three composite channels, in display order R, G, B."""
    return [spec for spec in CHANNELS if not spec.is_feature]


def feature_channels() -> list[ChannelSpec]:
    """The three narrowband panels, which are never blended into the RGB."""
    return [spec for spec in CHANNELS if spec.is_feature]


def mosaic_paths(
    outdir: str | Path, channels: list[str] | None = None
) -> dict[str, str]:
    """Map channel name to the mosaic path it should be read from."""
    table = channel_table(outdir)
    wanted = list(table) if channels is None else list(channels)
    missing = [name for name in wanted if name not in table]
    if missing:
        raise KeyError(f"unknown colour channel(s): {', '.join(missing)}")
    return {name: str(table[name].mosaic_path(outdir)) for name in wanted}


# --------------------------------------------------------------------------
# cross-filter registration
# --------------------------------------------------------------------------


@dataclass
class ChannelRegistration:
    """One channel's translation relative to the registration reference."""

    channel: str
    role: str
    is_reference: bool
    dx: float
    dy: float
    n_detected: int
    n_matched: int
    residual_rms_px: float
    residual_max_px: float
    offset_mad_px: float
    within_half_px: float
    median_arcsec: float
    negligible: bool

    @property
    def offset_px(self) -> float:
        return float(np.hypot(self.dx, self.dy))

    def as_dict(self) -> dict:
        return {
            "channel": self.channel,
            "role": self.role,
            "is_reference": self.is_reference,
            "dx_px": self.dx,
            "dy_px": self.dy,
            "offset_px": self.offset_px,
            "offset_arcsec": self.median_arcsec,
            "offset_mad_px": self.offset_mad_px,
            "residual_rms_px": self.residual_rms_px,
            "residual_max_px": self.residual_max_px,
            "n_detected": self.n_detected,
            "n_matched": self.n_matched,
            "within_0p5_px": self.within_half_px,
            "shift_negligible": self.negligible,
            "shift_applied": False,
        }


def _median_or_nan(values: np.ndarray) -> float:
    return float(np.median(values)) if np.size(values) else float("nan")


def solve_cross_filter_registration(
    catalogs: Mapping[str, tuple[np.ndarray, np.ndarray]],
    *,
    reference: str = COLOR_REFERENCE,
    match_radius_px: float = DEFAULT_MATCH_RADIUS_PX,
    grid_scale_arcsec: float = float("nan"),
    negligible_px: float = DEFAULT_NEGLIGIBLE_SHIFT_PX,
) -> list[ChannelRegistration]:
    """Match every channel to *reference* and solve a translation per channel.

    *catalogs* maps channel name to ``(x, y)`` in colour-grid pixels, as
    returned by :func:`jwst_stack.starcat.detect_tiled`.  The shift is the
    robust median of the per-star offsets (channel minus reference), the same
    estimator registration uses, so the number is comparable with the per-frame
    shifts in ``registration.json``.

    ``shift_applied`` is recorded as ``False`` for every channel and
    ``negligible`` is the real decision: a median below *negligible_px* is
    reported and left alone.  Correcting a sub-0.05 px difference measured
    against these same grids would be fitting the comparison floor, not fixing
    an error.
    """
    from jwst_stack.register import match_stars, solve_translation

    if reference not in catalogs:
        raise KeyError(f"registration reference {reference!r} is not in the catalogues")
    ref_x, ref_y = catalogs[reference]
    out: list[ChannelRegistration] = []
    for name, (x, y) in catalogs.items():
        spec = _CHANNEL_BY_NAME.get(name)
        role = spec.role if spec else "?"
        if name == reference:
            out.append(
                ChannelRegistration(
                    channel=name,
                    role=role,
                    is_reference=True,
                    dx=0.0,
                    dy=0.0,
                    n_detected=int(np.size(x)),
                    n_matched=int(np.size(x)),
                    residual_rms_px=float("nan"),
                    residual_max_px=float("nan"),
                    offset_mad_px=float("nan"),
                    within_half_px=float("nan"),
                    median_arcsec=0.0,
                    negligible=True,
                )
            )
            continue
        dx, dy, residuals = match_stars(x, y, ref_x, ref_y, match_radius_px)
        sx, sy, rms, mx = solve_translation(dx, dy, len(dx))
        offset = float(np.hypot(sx, sy)) if np.isfinite(sx) else float("nan")
        arcsec = (
            px_to_arcsec(offset, grid_scale_arcsec, GRID_PX)
            if np.isfinite(offset) and np.isfinite(grid_scale_arcsec)
            else float("nan")
        )
        out.append(
            ChannelRegistration(
                channel=name,
                role=role,
                is_reference=False,
                dx=float(sx),
                dy=float(sy),
                n_detected=int(np.size(x)),
                n_matched=int(len(dx)),
                residual_rms_px=float(rms),
                residual_max_px=float(mx),
                offset_mad_px=_median_or_nan(residuals),
                within_half_px=float(np.mean(residuals <= 0.5))
                if residuals.size
                else float("nan"),
                median_arcsec=float(arcsec),
                negligible=bool(np.isfinite(offset) and offset < negligible_px),
            )
        )
    return out


def measure_cross_filter_registration(
    outdir: str | Path,
    *,
    channels: list[str] | None = None,
    reference: str = COLOR_REFERENCE,
    tile_px: int = DEFAULT_REGISTER_TILE_PX,
    halo_px: int = DEFAULT_REGISTER_HALO_PX,
    match_radius_px: float = DEFAULT_MATCH_RADIUS_PX,
    negligible_px: float = DEFAULT_NEGLIGIBLE_SHIFT_PX,
    fwhm_px: float = DEFAULT_REGISTER_DETECT_FWHM_PX,
    progress=None,
) -> tuple[list[ChannelRegistration], float]:
    """Detect on every channel mosaic, then solve the cross-filter match.

    Streams each mosaic in tiles, so the 2.11 GB science array is never loaded
    whole.  Returns the per-channel solutions and the grid scale they were
    solved on.
    """
    from jwst_stack.starcat import covered_tiles, detect_tiled, open_mosaic

    paths = mosaic_paths(outdir, channels)
    if reference not in paths:
        raise KeyError(f"reference channel {reference!r} is not among {sorted(paths)}")
    catalogs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    scale = float("nan")
    for name, path in paths.items():
        data = open_mosaic(path)
        scale = data.grid_scale_arcsec
        tiles = covered_tiles(data.shape, data.coverage, tile_px)
        if progress is not None:
            progress(name, len(tiles))
        x, y, _ = detect_tiled(
            data.sci,
            tiles,
            halo_px=halo_px,
            fwhm_px=fwhm_px,
            max_stars=60000,
        )
        catalogs[name] = (x, y)
    solutions = solve_cross_filter_registration(
        catalogs,
        reference=reference,
        match_radius_px=match_radius_px,
        grid_scale_arcsec=scale,
        negligible_px=negligible_px,
    )
    return solutions, scale


def write_color_registration_json(
    path: str | Path,
    solutions: list[ChannelRegistration],
    *,
    reference: str = COLOR_REFERENCE,
    grid_scale_arcsec: float = float("nan"),
    grid_path: str | Path | None = None,
    match_radius_px: float = DEFAULT_MATCH_RADIUS_PX,
    negligible_px: float = DEFAULT_NEGLIGIBLE_SHIFT_PX,
    grids: dict[str, tuple[int, int]] | None = None,
) -> Path:
    """Write ``color_registration_f200w.json``.

    The offsets are colour-grid pixels, and the file says so, because a
    *cross-filter* offset has no single i2d scale to convert it: the channels
    are different filters and the measurement is in the grid they share.
    """
    from jwst_stack.io import update_json_section

    reference_solution = next((s for s in solutions if s.is_reference), None)
    others = [s for s in solutions if not s.is_reference]
    significant = [s for s in others if not s.negligible]
    payload = {
        "reference": reference,
        "pixel_basis": GRID_PX,
        "pixel_basis_label": px_label(GRID_PX),
        "grid_path": str(grid_path) if grid_path else None,
        "grid_scale_arcsec_per_px": float(grid_scale_arcsec),
        "match_radius_px": float(match_radius_px),
        "negligible_shift_px": float(negligible_px),
        "correction_applied": bool(significant),
        "correction_reason": (
            "one or more channels exceed the negligible-shift threshold and "
            "would need shifting before compositing"
            if significant
            else (
                "every channel's median shift is below the threshold, so the "
                "composite is built without any cross-filter shift: correcting "
                "a sub-threshold offset would be fitting the comparison floor"
            )
        ),
        # The reference channel's own angular pixel scale.  Every channel here
        # is resampled onto the same grid, so this is the grid scale; taking it
        # from the reference's *median shift* would be wrong, and silently so,
        # because the reference's shift is zero by definition.
        "reference_scale_arcsec_per_px": float(grid_scale_arcsec),
        "channels": {s.channel: s.as_dict() for s in solutions},
    }
    if grids:
        payload["grids"] = {k: list(v) for k, v in grids.items()}
    return update_json_section(path, "cross_filter_registration", payload)


def format_registration_summary_color(
    solutions: list[ChannelRegistration], grid_scale_arcsec: float = float("nan")
) -> str:
    """Terminal table for the cross-filter match."""
    lines = [
        f"cross-filter registration (reference = {COLOR_REFERENCE}), basis {px_label(GRID_PX)}"
    ]
    if np.isfinite(grid_scale_arcsec):
        lines.append(f"  grid scale: {grid_scale_arcsec:.4f} arcsec/px")
    lines.append(
        f"  {'channel':<14} {'role':<8} {'dx px':>9} {'dy px':>9} {'|t| px':>8} "
        f"{'|t| arcsec':>11} {'rms px':>7} {'matched':>9} {'<=0.5px':>8}"
    )
    for s in solutions:
        lines.append(
            f"  {s.channel:<14} {s.role:<8} {s.dx:+9.4f} {s.dy:+9.4f} "
            f"{s.offset_px:8.4f} {s.median_arcsec:11.5f} {s.residual_rms_px:7.3f} "
            f"{s.n_matched:9d} {s.within_half_px:8.3f}"
            if not s.is_reference
            else f"  {s.channel:<14} {s.role:<8} {'reference':>28}"
            f" {s.n_detected:9d}"
        )
    return "\n".join(lines)
