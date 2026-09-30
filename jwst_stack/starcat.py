"""Tiled star catalogues on a mosaic grid.

The measurements a colour composite needs before it can be built - how wide
each channel's PSF actually is, and how far each channel sits from the
registration reference - are both star measurements on the six finished
mosaics, and both have the same problem: the mosaics are 351.8 Mpx, so they
cannot be handed to a star finder whole.

This module holds the shared solution.  It streams the mosaic in tiles with a
halo, runs :func:`jwst_stack.register.detect_stars` (the same DAOStarFinder
settings registration and validation use, so the numbers are comparable) on
each halo-expanded tile, and keeps only the detections whose centroid falls
inside the tile's interior.  That interior test is what makes the result a
function of the mosaic rather than of the tiling: a star in the halo band of
one tile is rejected there and accepted in the next tile that owns it, so
every interior star is counted exactly once and the catalogue is independent
of ``tile_px``.

Two things about how a mosaic is stored are load-bearing and are centralised
here rather than re-derived at each call site, because getting either one
wrong is silent (see "FITS structure gotcha" in the project notes):

* the science array lives in the **primary** HDU, not in an ``SCI`` extension;
* the ``COVERAGE`` extension is an **unsigned** int16 depth - ``BITPIX 16`` with
  ``BZERO 32768`` - so it cannot be memory-mapped at all, and astropy refuses
  with ``ValueError: Cannot load a memory-mapped image``.  Its physical value
  is ``raw + 32768``, so even if it were mapped, ``> 0`` on the raw array would
  be true everywhere.  Coverage is therefore taken from the science array's NaN
  mask, which means the same thing (a pixel outside every footprint is NaN) and
  needs no second array: one 1.41 GB memmap instead of two.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from jwst_stack.mosaic import iter_tiles
from jwst_stack.register import detect_stars

#: Detection fwhm used when the caller has no better estimate.  A wrong value
#: only affects DAOStarFinder's detection efficiency, not the fitted widths in
#: :mod:`jwst_stack.psf`, so this is a starting point rather than a claim.
DEFAULT_DETECT_FWHM_PX = 3.0

#: Detection threshold in units of the robust (MAD) noise of the tile.
DEFAULT_THRESHOLD_SIGMA = 5.0

#: Halo added around each tile so a star near a tile edge is seen whole.
DEFAULT_HALO_PX = 24

#: Fraction of a tile that must be covered before it is worth detecting in.
DEFAULT_MIN_COVERAGE = 0.9

#: Stride used when sampling coverage, so tile selection never reads the whole
#: 351.8 Mpx array.
COVERAGE_STRIDE = 8

#: A query giving the covered fraction of one tile, ``probe(y0, y1, x0, x1)``.
CoverageProbe = Callable[[int, int, int, int], float]


@dataclass
class MosaicData:
    """A mosaic opened for measurement: science array, coverage query and WCS.

    ``sci`` is a memmapped view, not a copy - nothing here loads the 1.41 GB
    science array, and callers stream it tile by tile.
    """

    path: Path
    header: fits.Header
    wcs: WCS
    sci: np.ndarray
    coverage: CoverageProbe
    shape: tuple[int, int]
    grid_scale_arcsec: float

    @property
    def covered_pixels(self) -> int:
        """Covered pixel count, counted on a strided sample.

        Exact only at ``stride=1``; the mosaic is 351.8 Mpx, so the default
        is a 64x sample.  It is a report, not an input to any measurement.
        """
        step = COVERAGE_STRIDE
        sample = np.asarray(self.sci[::step, ::step])
        return int(np.count_nonzero(np.isfinite(sample)))


def _fraction(flags: np.ndarray) -> float:
    return float(np.count_nonzero(flags)) / float(flags.size) if flags.size else 0.0


def coverage_probe(
    sci: np.ndarray | None = None, coverage: np.ndarray | None = None
) -> CoverageProbe | None:
    """A covered-fraction query, or None when coverage is unknown.

    *coverage* is a depth array (``> 0`` is covered); *sci* is a science array
    whose NaNs are the holes.  Both are sampled on a stride, so the whole
    351.8 Mpx array is never read to choose tiles.
    """
    if coverage is not None:
        return lambda y0, y1, x0, x1: _fraction(
            np.asarray(coverage[y0:y1:COVERAGE_STRIDE, x0:x1:COVERAGE_STRIDE]) > 0
        )
    if sci is not None:
        return lambda y0, y1, x0, x1: _fraction(
            np.isfinite(sci[y0:y1:COVERAGE_STRIDE, x0:x1:COVERAGE_STRIDE])
        )
    return None


def open_mosaic(path: str | Path) -> MosaicData:
    """Open a mosaic for streaming measurement.

    The science array is memmapped and the coverage query is derived from its
    NaN mask, so this never materialises either array.  Callers should keep the
    returned object alive for as long as they read from it.
    """
    file = Path(path)
    hdul = fits.open(file, memmap=True)
    try:
        header = hdul[0].header
        sci = hdul[0].data
        if sci is None:
            raise ValueError(f"{file} has no science data in the primary HDU")
        wcs = WCS(header)
        shape = (int(sci.shape[0]), int(sci.shape[1]))
        probe = coverage_probe(sci=sci)
        if probe is None:  # pragma: no cover - sci is never None here
            raise ValueError(f"{file} has no science data to measure coverage from")
        scale = float(wcs.proj_plane_pixel_scales()[0].to_value("arcsec"))
    except BaseException:
        hdul.close()
        raise
    return MosaicData(
        path=file,
        header=header,
        wcs=wcs,
        sci=sci,
        coverage=probe,
        shape=shape,
        grid_scale_arcsec=scale,
    )


def tile_coverage_fraction(
    coverage: np.ndarray, y0: int, y1: int, x0: int, x1: int, stride: int = COVERAGE_STRIDE
) -> float:
    """Fraction of sampled pixels in the tile that are covered (depth > 0)."""
    return _fraction(np.asarray(coverage[y0:y1:stride, x0:x1:stride]) > 0)


def covered_tiles(
    shape: tuple[int, int],
    probe: CoverageProbe | None,
    tile_px: int,
    *,
    n_tiles: int | None = None,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
) -> list[tuple[int, int, int, int]]:
    """Interior tiles worth running a star finder on, in row-major order.

    With a *probe*, only tiles at least *min_coverage* covered are returned; a
    star catalogue is meaningless on a tile that is half hole.  Without one,
    every tile is returned, which is only sensible for a small array in a test.

    *n_tiles* caps the result by taking an evenly spaced subset of the covered
    tiles.  Sampling spreads the measurement over the field rather than
    concentrating it in one detector's corner, and bounds the runtime; pass
    ``None`` for a full-footprint catalogue.
    """
    tiles = iter_tiles(shape, tile_px)
    kept = list(tiles) if probe is None else [t for t in tiles if probe(*t) >= min_coverage]
    if n_tiles is None or n_tiles >= len(kept):
        return kept
    if n_tiles <= 0:
        return []
    # Evenly spaced indices over the row-major covered list, so the sample runs
    # across the field instead of clustering in the first few rows.
    picks = np.unique(np.linspace(0, len(kept) - 1, n_tiles).round().astype(int))
    return [kept[i] for i in picks]


def detect_tiled(
    image: np.ndarray,
    tiles: list[tuple[int, int, int, int]],
    *,
    halo_px: int = DEFAULT_HALO_PX,
    keep_margin: float = 0.0,
    fwhm_px: float = DEFAULT_DETECT_FWHM_PX,
    threshold_sigma: float = DEFAULT_THRESHOLD_SIGMA,
    per_tile_stars: int = 400,
    max_stars: int | None = None,
    progress=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Detect stars over *tiles*, returned in whole-image pixel coordinates.

    Each tile is expanded by *halo_px* before detection, so DAOStarFinder sees
    a star whole even when its centroid is near the tile boundary.  Only
    centroids inside the un-expanded interior - inset by a further
    *keep_margin* - are kept, so:

    * no star is counted twice, whatever the halo, because the overlap band
      belongs to the neighbouring tile's interior and not to this one;
    * a caller that needs room to cut a patch around each centroid (the FWHM
      fit) can ask for it with *keep_margin* and still get every star.

    *per_tile_stars* bounds each individual DAOStarFinder call; *max_stars*
    bounds the returned catalogue as a whole, and defaults to no bound.  They
    are separate because a full-footprint run wants thousands of stars while a
    single-tile brightness ranking does not.

    Returns ``(x, y, flux)`` sorted by decreasing flux.
    """
    arr = np.asarray(image)
    ny, nx = arr.shape
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    fs: list[np.ndarray] = []
    for i, (y0, y1, x0, x1) in enumerate(tiles, start=1):
        ey0 = max(0, y0 - halo_px)
        ex0 = max(0, x0 - halo_px)
        ey1 = min(ny, y1 + halo_px)
        ex1 = min(nx, x1 + halo_px)
        patch = np.ascontiguousarray(arr[ey0:ey1, ex0:ex1], dtype=np.float32)
        dx, dy, df = detect_stars(
            patch,
            fwhm_px=fwhm_px,
            threshold_sigma=threshold_sigma,
            max_stars=per_tile_stars,
        )
        if dx.size:
            keep = (
                (dx >= (x0 - ex0) + keep_margin)
                & (dx < (x1 - ex0) - keep_margin)
                & (dy >= (y0 - ey0) + keep_margin)
                & (dy < (y1 - ey0) - keep_margin)
            )
            if np.any(keep):
                xs.append(dx[keep] + ex0)
                ys.append(dy[keep] + ey0)
                fs.append(df[keep])
        if progress is not None:
            progress(i, len(tiles))
    if not xs:
        return (
            np.empty(0, dtype=float),
            np.empty(0, dtype=float),
            np.empty(0, dtype=float),
        )
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    f = np.concatenate(fs)
    order = np.argsort(f)[::-1]
    if max_stars is not None:
        order = order[:max_stars]
    return x[order], y[order], f[order]


def brightest_separated(
    x: np.ndarray, y: np.ndarray, flux: np.ndarray, count: int, separation: float
) -> np.ndarray:
    """Indices of the *count* brightest stars at least *separation* apart.

    Greedy on a brightness-sorted list, rejecting anything within
    *separation* of an already-accepted star.  This is how the FWHM sample is
    kept independent: a blend or a saturated double core is rejected as a
    neighbour rather than averaged into the channel's PSF width.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size == 0 or count <= 0 or separation <= 0:
        return np.empty(0, dtype=int)
    order = np.argsort(np.asarray(flux, dtype=float))[::-1]
    accepted: list[int] = []
    ax: list[float] = []
    ay: list[float] = []
    for idx in order:
        px = float(x[idx])
        py = float(y[idx])
        if any((px - qx) ** 2 + (py - qy) ** 2 < separation**2 for qx, qy in zip(ax, ay)):
            continue
        accepted.append(int(idx))
        ax.append(px)
        ay.append(py)
        if len(accepted) >= count:
            break
    return np.asarray(accepted, dtype=int)
