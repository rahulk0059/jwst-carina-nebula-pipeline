"""Phase 4: the tiled colour builder.

Turns the six co-located mosaics into the shipped colour products - one RGB
cube (F444W;CLEAR / F200W / F090W) and three single-band feature panels (F187N,
F335M, F444W;F470N) - by streaming output tiles and never holding a full
channel in memory.

The two honesty decisions live in :mod:`jwst_stack.color` and are applied here
unchanged: each channel is background-subtracted against its *own*
:func:`~jwst_stack.color.source_free_level`, and all six then share one
``lo``/``hi``/softening triple from :func:`~jwst_stack.color.shared_limits`.  A
per-channel min-max would rescale a faint channel onto a bright one and destroy
the colour; a shared stretch keeps relative brightness, which is the whole point.

Why the builder is tiled rather than one big call to
:func:`jwst_stack.color.stretch_channels`.  The grid is 15895 x 22130 = 351.8
Mpx, so a float32 array per channel is 1.41 GB and the three RGB bands alone are
4.22 GB.  :func:`stretch_channels` is the readable statement of the policy and
the target of the invariance tests; it is not something to point at 2.11 GB
files.  Here the statistics are measured once over the memmapped mosaics
(:func:`plan_color`) and then applied tile by tile (:func:`build_color`).

Pixel basis is ``grid_px`` throughout - every mosaic is on the shared
0.031 arcsec grid, so there is no resampling and no second scale to get wrong.
The cross-filter shifts measured in Phase 2 are *not* applied: all five medians
came out below the 0.25 px threshold, so ``correction_applied`` is false and
shifting by a sub-threshold offset would be fitting the comparison floor.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits

from jwst_stack import color
from jwst_stack.color import (
    CHANNELS,
    DEFAULT_CLIP_SIGMA,
    DEFAULT_HI_PCT,
    DEFAULT_LO_PCT,
    DEFAULT_SOFTENING_FRACTION,
    DEFAULT_TILE_PX,
    ChannelSpec,
    background_centred_asinh,
    default_softening,
    shared_limits,
    source_free_level,
    subtract_background,
)
from jwst_stack.grid import load_grid
from jwst_stack.mosaic import iter_tiles, peak_rss_gb
from jwst_stack.starcat import open_mosaic

#: Rows written per streaming step.  Full-width, so a band of a C-ordered array
#: is contiguous and can be handed to the FITS writer verbatim - see
#: :func:`_stream_products` for why the decomposition is rows and not 2-D tiles.
#: 128 rows x 22130 cols is 11.3 MB per channel, so six channels plus the uint16
#: outputs is a ~100 MB working set regardless of the grid size.
DEFAULT_ROW_PX = 128

#: The uint16 ceiling the products are written at.  The composite is a display
#: product, so 16 bits is ample and halves the output against 32-bit float.
UINT16_MAX = 65535

#: Longest edge of the PNG preview, in pixels.  A full-plane PNG would be ~1 GB
#: and unreadable; the FITS is the product, the PNG is for looking at.
DEFAULT_PREVIEW_PX = 2000

#: Pixels sampled per channel when pooling the shared limits.  Each mosaic is
#: 351.8 Mpx; half a million evenly spaced finite samples pin a percentile to
#: well inside its own sampling error and cost 2 MB per channel.
POOL_SAMPLE_PX = 500_000


@dataclass
class ColorPlan:
    """Everything needed to write the products, measured from the mosaics.

    Produced by :func:`plan_color` and reusable across several builds - it is
    the expensive half, and rebuilding it with a different tile size must not
    change the answer, which the tests check.
    """

    grid_path: Path
    grid_scale_arcsec: float
    shape: tuple[int, int]
    backgrounds: dict[str, float]
    lo: float
    hi: float
    softening: float
    lo_pct: float = DEFAULT_LO_PCT
    hi_pct: float = DEFAULT_HI_PCT
    softening_fraction: float = DEFAULT_SOFTENING_FRACTION
    channels: tuple[ChannelSpec, ...] = CHANNELS

    def as_dict(self) -> dict:
        return {
            "pixel_basis": "grid_px",
            "grid_path": str(self.grid_path),
            "grid_scale_arcsec_per_px": self.grid_scale_arcsec,
            "shape": [int(self.shape[0]), int(self.shape[1])],
            "background_policy": "self-referenced per channel (additive, sigma-clipped median)",
            "backgrounds": {k: float(v) for k, v in self.backgrounds.items()},
            "stretch_policy": (
                "one shared lo/hi/softening pooled over all six channels after "
                "per-channel background subtraction"
            ),
            "lo": self.lo,
            "hi": self.hi,
            "softening": self.softening,
            "lo_pct": self.lo_pct,
            "hi_pct": self.hi_pct,
            "softening_fraction": self.softening_fraction,
            "note": (
                "a shared stretch preserves relative brightness between channels; "
                "it does not equalise them, so a channel with a larger residual "
                "nebular gradient occupies more of the display range"
            ),
        }


@dataclass
class ColorBuildResult:
    """What a build produced, and what it cost."""

    plan: ColorPlan
    products: list[Path] = field(default_factory=list)
    previews: list[Path] = field(default_factory=list)
    n_tiles: int = 0
    covered_pixels: int = 0
    wall_seconds: float = 0.0
    peak_rss_gb: float = 0.0
    tile_px: int = DEFAULT_ROW_PX

    def as_dict(self) -> dict:
        return {
            "products": [str(p) for p in self.products],
            "previews": [str(p) for p in self.previews],
            "n_tiles": self.n_tiles,
            "covered_pixels": self.covered_pixels,
            "wall_seconds": round(self.wall_seconds, 2),
            "peak_rss_gb": round(self.peak_rss_gb, 3),
            "tile_px": self.tile_px,
            "byte_size_estimate": self.estimate_bytes(),
        }

    def estimate_bytes(self) -> dict[str, int]:
        """Predicted on-disk size of each product, from shape and dtype.

        Exact, not a guess: a uint16 product's size is fully determined by its
        shape, so this is what the disk check before a full build should compare
        against.
        """
        ny, nx = self.plan.shape
        sizes = {
            "rgb_cube_uint16": 3 * ny * nx * 2,
            "each_feature_uint16": ny * nx * 2,
            "all_four_uint16": (3 + 3) * ny * nx * 2,
        }
        return sizes


def _open(specs, outdir):
    return {spec.name: open_mosaic(spec.mosaic_path(outdir)) for spec in specs}


def plan_color(
    outdir: str | Path,
    grid: str | Path,
    channels: list[str] | None = None,
    *,
    tile_px: int = DEFAULT_TILE_PX,
    lo_pct: float = DEFAULT_LO_PCT,
    hi_pct: float = DEFAULT_HI_PCT,
    softening_fraction: float = DEFAULT_SOFTENING_FRACTION,
    clip_sigma: float = DEFAULT_CLIP_SIGMA,
) -> ColorPlan:
    """Measure the per-channel backgrounds and the shared stretch from the mosaics.

    Reads all six mosaics twice: once for each channel's robust level, once
    strided for the pooled limits.  Both passes are memmap-friendly by
    construction - :func:`source_free_level` walks tiles taking local
    statistics, and the limits are pooled from a strided sample - so neither
    materialises a channel.

    No i2d is consulted and no channel is scaled against another.
    """
    outdir = Path(outdir)
    table = color.channel_table(outdir)
    names = list(channels) if channels else [spec.name for spec in CHANNELS]
    specs = [table[n] for n in names]

    wcs, shape = load_grid(Path(grid))
    scale = float(wcs.proj_plane_pixel_scales()[0].to_value("arcsec"))

    backgrounds: dict[str, float] = {}
    samples: list[np.ndarray] = []
    for spec in specs:
        data = open_mosaic(spec.mosaic_path(outdir))
        try:
            backgrounds[spec.name] = source_free_level(
                data.sci, tile_px=tile_px, clip_sigma=clip_sigma
            )
            flat = np.asarray(data.sci).ravel()
            step = max(1, flat.size // POOL_SAMPLE_PX)
            sample = np.asarray(flat[::step], dtype=float)
            sample = sample[np.isfinite(sample)] - backgrounds[spec.name]
            if sample.size:
                samples.append(sample)
        finally:
            # Dropping the handle closes the memmap.  Six open 1.41 GB memmaps
            # at once is page-cache pressure for no benefit: the statistics are
            # taken here, not streamed.
            del data
    lo, hi = shared_limits(samples, lo_pct=lo_pct, hi_pct=hi_pct)
    return ColorPlan(
        grid_path=Path(grid),
        grid_scale_arcsec=scale,
        shape=(int(shape[0]), int(shape[1])),
        backgrounds=backgrounds,
        lo=lo,
        hi=hi,
        softening=default_softening(lo, hi, softening_fraction),
        lo_pct=lo_pct,
        hi_pct=hi_pct,
        softening_fraction=softening_fraction,
        channels=tuple(specs),
    )


def _stretch(tile: np.ndarray, plan: ColorPlan, channel: str) -> np.ndarray:
    """Subtract *channel*'s own background, then asinh-stretch to [0, 1].

    The subtraction is per channel and from the plan, never a shared constant
    and never anything read from an i2d - that is the whole content of the
    self-referenced-background policy, and it is applied here at the last
    possible moment so it cannot be forgotten by a caller.

    Delegates to the same :func:`~jwst_stack.color.background_centred_asinh` the
    policy tests exercise, with ``background=0.0`` because the pedestal has
    already been removed - so a row band inside a build is stretched by exactly
    the function a whole array uses in a test.
    """
    return background_centred_asinh(
        subtract_background(tile, plan.backgrounds[channel]),
        background=0.0,
        a=plan.softening,
        lo=plan.lo,
        hi=plan.hi,
    )


def _to_uint16(stretched: np.ndarray) -> np.ndarray:
    """Scale [0, 1] to uint16, writing uncovered pixels as 0.

    0 is the right fill for an uncovered pixel: the stretch maps a
    background-subtracted background to somewhere near the middle of the range,
    so 0 is distinguishable from every real value, and a NaN cannot be stored.
    """
    out = np.zeros(stretched.shape, dtype=np.uint16)
    finite = np.isfinite(stretched)
    out[finite] = np.rint(
        np.clip(stretched[finite], 0.0, 1.0) * UINT16_MAX
    ).astype(np.uint16)
    return out


#: FITS ``BZERO`` for the uint16 products, set by astropy's own convention.
FITS_BZERO_UINT16 = 32768


def _to_stored_int16(values: np.ndarray) -> np.ndarray:
    """uint16 -> the int16 bytes a FITS ``BZERO 32768`` extension actually stores.

    FITS defines ``physical = BZERO + BSCALE * stored``, so the stored value is
    ``physical - 32768``.  That is an arithmetic shift, **not** a bit
    reinterpretation: ``values.view(np.int16)`` looks like it should work and
    does not, and it fails silently - it writes 32768 where 0 was meant and
    comes back as an image full of 32768.  Every streamed chunk goes through
    here.
    """
    return (values.astype(np.int32) - FITS_BZERO_UINT16).astype(np.int16)


def _uint16_header(shape: tuple[int, ...], plan: ColorPlan, bandpasses: list[str]) -> fits.Header:
    """A FITS header for a uint16 product, built from astropy's own convention.

    FITS has no unsigned 16-bit type, so a uint16 array is stored as ``BITPIX
    16`` with ``BZERO 32768``.  astropy encodes that convention on a throwaway
    HDU here rather than hand-rolling it, and the data is converted by
    :func:`_to_stored_int16`; this is the same trap as the mosaics' ``COVERAGE``
    extension, which is why coverage is read from the science NaN mask instead
    of from that array.
    """
    prototype = fits.PrimaryHDU(data=np.zeros(shape, dtype=np.uint16)).header
    header = prototype.copy()
    header["BUNIT"] = "DN/s/px"
    header["PIXELAT"] = "grid_px"
    header["PIXSCALE"] = (plan.grid_scale_arcsec, "arcsec/pixel")
    for index, bandpass in enumerate(bandpasses, start=1):
        header[f"BAND{index}"] = (bandpass, f"band {index} filter")
    return header


def _stream_products(outdir, plan, data, rgb_path, feature_path, row_px, progress):
    """Write the cube and the panels a full-width row band at a time.

    Why rows and not 2-D tiles.  ``fits.StreamingHDU`` writes the buffer it is
    handed *verbatim*, so every chunk must be C-contiguous.  In a C-ordered
    array only a full-width row range is contiguous - an arbitrary
    ``[y0:y1, x0:x1]`` block is not - so the decomposition is rows and a 2-D
    tile would need ``ascontiguousarray`` on every channel, every row.  The
    cost is bounded memory: one row band of six float32 channels is ~11 MB at
    128 rows on this grid.

    The products deliberately use different in-file orders because only
    ``(ny, nx)`` row bands stream correctly: the feature panels are genuine
    ``(ny, nx)`` images, and the RGB cube is ``(ny, nx, 3)`` - band **last** -
    so every row band is a contiguous ``(row_px, nx, 3)`` block and a write is
    still one buffer.  astropy reads both back as the natural array; in
    particular the cube is *not* band-major ``(3, ny, nx)``, because a
    partial 3-D write through ``fits.StreamingHDU`` does not land where it
    claims (verified empirically: a ``(3, y, nx)`` band is written misaligned,
    a ``(y, nx, 3)`` one is exact).  The band order is unambiguously recorded
    in the ``BAND1/2/3`` header cards.

    **The cube masks its holes across channels, not per channel.**  Each RGB
    band individually writes 0 where *its own* mosaic has no data, so a pixel
    left uncovered by one footprint edge (for example F200W's, which the real
    field shows has a partial fringe against F090W/F444W) would survive as
    ``(R, 0, B)`` - bright magenta, two channels' light with the third's
    footprint cut out.  A pixel is therefore zeroed in **all three** bands
    whenever any one of the RGB mosaics lacks data there, so the composite
    never shows light from a hole.  Coverage is read from the raw band - a
    value that is finite and non-zero - not from the stretched result, where a
    naked 0 is indistinguishable from an uncovered pixel.  The feature panels
    are single-band below and keep the per-channel fill.
    """
    ny, nx = plan.shape
    table = color.channel_table(outdir)
    rgb_names = [s.name for s in plan.channels if not s.is_feature]
    feature_names = [s.name for s in plan.channels if s.is_feature]

    cube_hdu = fits.StreamingHDU(
        rgb_path,
        header=_uint16_header((ny, nx, 3), plan, [table[n].bandpass for n in rgb_names]),
    )
    feat_hdus = {
        name: fits.StreamingHDU(
            feature_path[name],
            header=_uint16_header((ny, nx), plan, [table[name].bandpass]),
        )
        for name in feature_names
    }
    covered = 0
    n_bands_written = 0
    rgb_index = {name: i for i, name in enumerate(rgb_names)}
    try:
        for y0 in range(0, ny, row_px):
            y1 = min(ny, y0 + row_px)
            plane = np.empty((y1 - y0, nx, 3), dtype=np.uint16)
            rgb_cov = np.ones((y1 - y0, nx), dtype=bool)
            for spec in plan.channels:
                chunk = np.asarray(data[spec.name].sci[y0:y1, :], dtype=float)
                cov = np.isfinite(chunk) & (chunk != 0)
                if spec.name == rgb_names[0]:
                    covered += int(np.count_nonzero(np.isfinite(chunk)))
                written = _to_uint16(_stretch(chunk, plan, spec.name))
                if spec.is_feature:
                    feat_hdus[spec.name].write(_to_stored_int16(written))
                else:
                    plane[:, :, rgb_index[spec.name]] = written
                    rgb_cov &= cov
            plane[~rgb_cov] = 0
            cube_hdu.write(_to_stored_int16(plane))
            n_bands_written += 1
            if progress and n_bands_written % 50 == 0:
                print(f"  {y1}/{ny} rows")
    finally:
        cube_hdu.close()
        for hdu in feat_hdus.values():
            hdu.close()
    return covered, n_bands_written


def build_color(
    outdir: str | Path,
    plan: ColorPlan,
    *,
    rgb_out: str | Path | None = None,
    feature_out: dict[str, str | Path] | None = None,
    row_px: int = DEFAULT_ROW_PX,
    previews: bool = True,
    preview_px: int = DEFAULT_PREVIEW_PX,
    progress: bool = False,
) -> ColorBuildResult:
    """Write the RGB cube and the feature panels by streaming full-width row bands.

    Each band is read from all six memmapped mosaics, stretched, and handed to
    the streaming writer, so the only full-size arrays in play are the six
    memmaps (page cache, not heap) and the writers' own buffers.

    *plan* is passed in rather than measured here so the statistics are computed
    once and every product in a run is stretched identically.
    """
    started = time.time()
    outdir = Path(outdir)
    table = color.channel_table(outdir)
    rgb_names = [s.name for s in plan.channels if not s.is_feature]
    feature_names = [s.name for s in plan.channels if s.is_feature]
    rgb_path = Path(rgb_out) if rgb_out else outdir / "color_rgb.fits"
    feature_path = {
        name: Path((feature_out or {}).get(name, outdir / f"color_{_stem(name)}.fits"))
        for name in feature_names
    }

    result = ColorBuildResult(plan=plan, tile_px=row_px)
    for path in [rgb_path, *feature_path.values()]:
        path.unlink(missing_ok=True)
    data = _open(plan.channels, outdir)
    try:
        covered, n_bands = _stream_products(
            outdir, plan, data, rgb_path, feature_path, row_px, progress
        )
    finally:
        data.clear()

    result.covered_pixels = covered
    result.n_tiles = n_bands
    result.products = [rgb_path, *[feature_path[n] for n in feature_names]]
    if previews:
        result.previews = [
            _write_preview(rgb_path, outdir / "color_rgb.png", preview_px, ("R", "G", "B")),
            *[
                _write_preview(
                    feature_path[name],
                    outdir / f"color_{_stem(name)}.png",
                    preview_px,
                    (table[name].bandpass,),
                )
                for name in feature_names
            ],
        ]
    result.wall_seconds = time.time() - started
    result.peak_rss_gb = peak_rss_gb()
    return result


def _stem(name: str) -> str:
    """``F444W;F470N`` -> ``f470n``; the pupil is the distinguishing part."""
    parts = [p for p in name.split(";") if p]
    return parts[-1].lower() if len(parts) > 1 else parts[0].lower()


def _write_preview(
    path: str | Path, out_png: Path, target: int = DEFAULT_PREVIEW_PX,
    labels: tuple[str, ...] = ("",),
) -> Path:
    """Block-reduce a written product to a PNG of at most *target* px per axis.

    The reduction is a nan-aware **mean** over integer blocks, so surface
    brightness is preserved rather than resampled, and no interpolation can mix
    neighbouring pixels - the point of the preview is to show what was written,
    not to look smooth.  Both products are handled: the ``(ny, nx)`` feature
    panels, rendered as a single intensity image, and the ``(ny, nx, 3)``
    band-last cube, rendered as **one merged true-colour image** (the three
    reduced bands stacked into an RGB array for ``imshow``, not three side by
    side panels); because the product is already stretched by the shared
    lo/hi/softening, normalising to 0-1 keeps the relative brightness honest.

    A scaled image cannot be memory-mapped (astropy refuses when BZERO is
    present), and reading it back *section by section* is the wrong answer too:
    the quantum is 2.11 GB and each ``fits.getdata(section=...)`` re-opens and
    seeks, which turns the preview into tens of thousands of opens.  The
    product is therefore loaded once - this is the one place a full product is
    materialised, on purpose, after the streaming build has released it - and
    the mean is computed vectorised over integer blocks.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with fits.open(path, memmap=False) as hdul:
        header = hdul[0].header
        scaled = np.asarray(hdul[0].data)
    band_last = scaled.ndim == 3 and scaled.shape[-1] == 3
    n_bands = 3 if band_last else 1
    ny, nx = scaled.shape[0], scaled.shape[1]
    fy = max(1, ny // target)
    fx = max(1, nx // target)
    rows = max(1, ny // fy)
    cols = max(1, nx // fx)
    small = np.full((n_bands, rows, cols), np.nan, dtype=np.float64)
    for r in range(rows):
        y0, y1 = r * fy, min(ny, (r + 1) * fy)
        for c in range(cols):
            x0, x1 = c * fx, min(nx, (c + 1) * fx)
            block = scaled[y0:y1, x0:x1].astype(np.float64)
            with np.errstate(invalid="ignore"):
                if band_last:
                    small[:, r, c] = np.nanmean(block, axis=(0, 1))
                else:
                    small[0, r, c] = np.nanmean(block)
    del scaled
    scaled = np.nan_to_num(small / UINT16_MAX, nan=0.0)
    if band_last:
        rgb = np.transpose(scaled, (1, 2, 0))
        bands = [header.get(f"BAND{i}") for i in (1, 2, 3)]
        title = (
            f"R={bands[0]}  G={bands[1]}  B={bands[2]}"
            if all(bands)
            else "R  G  B"
        )
        fig, ax = plt.subplots(1, 1, figsize=(4.8, 4.4))
        ax.imshow(rgb, origin="lower", interpolation="nearest")
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        fig, axes = plt.subplots(1, n_bands, figsize=(4.2 * n_bands, 4.4), squeeze=False)
        for b in range(n_bands):
            axes[0][b].imshow(
                scaled[b], origin="lower", interpolation="nearest", vmin=0.0, vmax=1.0
            )
            axes[0][b].set_title(labels[b] if b < len(labels) else "", fontsize=9)
            axes[0][b].set_xticks([])
            axes[0][b].set_yticks([])
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return out_png


def export_rgb_viewer(
    rgb_path: str | Path,
    viewer_path: str | Path | None = None,
    progress: bool = False,
) -> Path:
    """Write the band-last cube as a standard band-first ``(3, ny, nx)`` cube.

    ``color_rgb.fits`` is intentionally band-last ``(ny, nx, 3)`` because that is
    the only decomposition ``fits.StreamingHDU`` lands row bands in exactly
    (a band-first chunk is written misaligned - see :func:`build_color`).  That
    is a *storage* choice; several viewers and editors - Siril being the one
    this export exists for - expect the conventional FITS RGB cube, where the
    colour axis is the slowest (``NAXIS1=nx, NAXIS2=ny, NAXIS3=3``) and each
    band occupies one contiguous plane.  This exports exactly that: the same
    values, band-major, in ``color_rgb_viewer.fits``.

    Nothing about the band-last file changes; this is a reader-oriented copy of
    its bytes.  The export streams one full-width band plane at a time so it
    never materialises the whole cube: the memmapped source is sliced per band
    and each ``(ny, nx)`` plane is written in one buffer (a plane is contiguous
    in the band-first layout, so the docstring's partial-3D-write caveat does
    not apply to it).  Headers are carried over where they mean something to a
    viewer - ``BUNIT``, ``PIXELAT``, ``PIXSCALE``, ``BAND1/2/3`` - plus
    ``CTYPE3='BAND'`` marking the colour axis.  The output is replaced (unlinked
    first), so a re-export cannot append to a stale file.

    Returns the ``viewer_path`` that was written.
    """
    rgb_path = Path(rgb_path)
    viewer_path = (
        Path(viewer_path)
        if viewer_path is not None
        else rgb_path.with_name("color_rgb_viewer.fits")
    )
    if viewer_path.exists():
        viewer_path.unlink()

    with fits.open(rgb_path, memmap=True, do_not_scale_image_data=True) as src:
        stored = src[0].data  # int16 (ny, nx, 3) memmap, no BZERO applied
        ny, nx, nband = stored.shape
        if nband != 3:
            raise ValueError(f"{rgb_path} has {nband} bands; this export is RGB-only")
        header = src[0].header
        bandpasses = [header.get(f"BAND{i}", "") for i in range(1, nband + 1)]

        prototype = fits.PrimaryHDU(data=np.zeros((nband, ny, nx), dtype=np.uint16)).header
        hdr = prototype.copy()
        for card in ("BUNIT", "PIXELAT", "PIXSCALE"):
            if card in header:
                hdr[card] = header[card]
        for i, bandpass in enumerate(bandpasses, start=1):
            hdr[f"BAND{i}"] = (bandpass, f"band {i} filter")
        hdr["CTYPE3"] = ("BAND", "RGB colour axis (band-first)")

        sh = fits.StreamingHDU(viewer_path, header=hdr)
        try:
            for b in range(nband):
                plane = np.ascontiguousarray(stored[:, :, b])
                sh.write(plane)
                if progress:
                    print(f"  band {b + 1} ({bandpasses[b]}): {ny}x{nx} plane")
        finally:
            sh.close()
    return viewer_path


def write_color_build_json(
    path: str | Path, result: ColorBuildResult, extra: dict | None = None
) -> Path:
    """Merge the build's plan and cost into the ``color_build`` JSON section."""
    from jwst_stack.io import update_json_section

    payload = {"plan": result.plan.as_dict(), "build": result.as_dict()}
    if extra:
        payload.update(extra)
    return update_json_section(path, "color_build", payload)
