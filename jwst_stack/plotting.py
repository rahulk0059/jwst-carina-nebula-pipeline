"""PNG previews and comparison of a stack against an official i2d mosaic."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from scipy.ndimage import maximum_filter
from scipy.stats import median_abs_deviation

from jwst_stack.align import reproject_to_grid


def asinh_stretch(image: np.ndarray, lo_pct: float = 1.0, hi_pct: float = 99.5) -> np.ndarray:
    """Return a normalized asinh-stretched image for display.

    Percentile clipping is applied first; the result is mapped to [0, 1].
    """
    finite = image[np.isfinite(image)]
    if finite.size == 0:
        return np.zeros_like(image, dtype=float)
    lo, hi = np.percentile(finite, [lo_pct, hi_pct])
    if hi <= lo:
        return np.zeros_like(image, dtype=float)
    center = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    stretched = np.arcsinh((image - center) / half)
    out = np.nan_to_num(stretched, nan=float(np.nanmin(stretched)))
    out -= np.nanmin(out)
    out /= np.nanmax(out)
    return out


def save_preview_png(
    image: np.ndarray,
    out_path: str | Path,
    title: str = "",
    dpi: int = 150,
) -> Path:
    """Write an asinh-stretched preview PNG of *image*."""
    out_path = Path(out_path)
    disp = asinh_stretch(image)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(disp, origin="lower", cmap="Greys_r", interpolation="nearest")
    if title:
        ax.set_title(title)
    ax.set_xlabel("output pixel")
    ax.set_ylabel("output pixel")
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


def detect_stars(
    image: np.ndarray, size: int = 11, n: int = 30
) -> np.ndarray:
    """Simple local-maximum star positions (pixel coords), brightest first."""
    finite = image[np.isfinite(image)]
    if finite.size == 0:
        return np.empty((0, 2))
    threshold = np.percentile(finite, 99.0)
    maxima = maximum_filter(image, size=size)
    is_peak = (image == maxima) & (image >= threshold)
    ys, xs = np.where(is_peak)
    flux = image[ys, xs]
    order = np.argsort(flux)[::-1][:n]
    return np.column_stack((xs[order], ys[order]))


def reproject_i2d(i2d_path: str | Path, output_wcs: WCS) -> np.ndarray:
    """Reproject an official ``*_i2d.fits`` mosaic SCI extension onto our grid."""
    with fits.open(i2d_path) as hdul:
        data = hdul["SCI"].data
        wcs = WCS(hdul["SCI"].header)
    on_grid, _ = reproject_to_grid(data, wcs, output_wcs, method="interp")
    return np.asarray(on_grid)


@dataclass
class ComparisonMetrics:
    """Simple agreement metrics between my stack and the official mosaic."""

    overlap_pixels: int
    median_diff: float
    mad_diff: float
    rms_diff: float
    matched_stars: int
    median_star_offset_px: float
    star_matches_a: np.ndarray
    star_matches_b: np.ndarray


def compare_stacks(
    mine: np.ndarray,
    official: np.ndarray,
    star_radius_px: float = 8.0,
) -> ComparisonMetrics:
    """Report background level and stellar agreement between two images."""
    both = np.isfinite(mine) & np.isfinite(official)
    diff = mine[both] - official[both]
    metrics = ComparisonMetrics(
        overlap_pixels=int(both.sum()),
        median_diff=float(np.median(diff)),
        mad_diff=float(median_abs_deviation(diff)),
        rms_diff=float(np.sqrt(np.mean(diff**2))),
        matched_stars=0,
        median_star_offset_px=float("nan"),
        star_matches_a=np.empty((0, 2)),
        star_matches_b=np.empty((0, 2)),
    )
    stars_mine = detect_stars(np.where(both, mine, np.nan))
    stars_official = detect_stars(np.where(both, official, np.nan))
    if len(stars_official):
        dist, idx = _nearest(stars_official, stars_mine)
        matched = dist <= star_radius_px
        if np.any(matched):
            metrics.star_matches_a = stars_official[matched]
            metrics.star_matches_b = stars_mine[idx[matched]]
            metrics.matched_stars = int(matched.sum())
            metrics.median_star_offset_px = float(np.median(dist[matched]))
    return metrics


def _nearest(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Distances and indices of the nearest point in *b* for each point in *a*."""
    if len(b) == 0:
        return np.array([np.inf] * len(a)), np.zeros(len(a), dtype=int)
    from scipy.spatial import cKDTree

    dist, idx = cKDTree(b).query(a, k=1)
    return dist, idx


def save_difference_png(
    mine: np.ndarray,
    official: np.ndarray,
    out_path: str | Path,
    title: str = "",
) -> Path:
    """Save a preview of their difference (NaN where only one image is valid)."""
    out_path = Path(out_path)
    diff = mine - official
    disp = asinh_stretch(np.ma.filled(np.ma.masked_invalid(diff), 0.0))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(disp, origin="lower", cmap="RdBu_r", interpolation="nearest", vmin=-1, vmax=1)
    if title:
        ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path