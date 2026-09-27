from __future__ import annotations

import numpy as np
from astropy.wcs import WCS


def make_wcs(
    crval: tuple[float, float] = (160.0, -58.0),
    crpix: tuple[float, float] = (32.5, 32.5),
    cd: np.ndarray | None = None,
) -> WCS:
    """A simple TAN WCS at NIRCam-like 0.031 arcsec/pixel."""
    if cd is None:
        cd = np.array([[-8.6e-6, 0.0], [0.0, 8.6e-6]])
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ["deg", "deg"]
    w.wcs.crpix = crpix
    w.wcs.crval = crval
    w.wcs.cd = cd
    return w


def gaussian_frame(
    wcs: WCS, peak_px: tuple[float, float], sigma: float = 1.5, amp: float = 100.0
) -> np.ndarray:
    """A 64x64 image with a Gaussian peak at ``peak_px`` (pixel, 0-indexed)."""
    ny, nx = 64, 64
    y, x = np.mgrid[0:ny, 0:nx]
    return amp * np.exp(
        -((x - peak_px[0]) ** 2 + (y - peak_px[1]) ** 2) / (2 * sigma**2)
    )


def shift_crval_to_match(wcs: WCS, peak_own: tuple[float, float], peak_target: tuple[float, float]) -> WCS:
    """Return a copy of *wcs* whose world position at ``peak_target`` equals
    the original world position at ``peak_own`` (same sky point, different
    pixel location).

    The shift is solved in world coordinates using the WCS's own transforms.
    A CD-based linear guess would be wrong for RA on a TAN projection by a
    factor of 1/cos(dec).
    """
    shifted = wcs.deepcopy()
    ra_own, dec_own = wcs.all_pix2world(peak_own[0], peak_own[1], 0)
    world_own = np.array([float(np.atleast_1d(ra_own)[0]), float(np.atleast_1d(dec_own)[0])])
    ra_tgt, dec_tgt = wcs.all_pix2world(peak_target[0], peak_target[1], 0)
    world_tgt = np.array([float(np.atleast_1d(ra_tgt)[0]), float(np.atleast_1d(dec_tgt)[0])])
    shifted.wcs.crval = np.array(shifted.wcs.crval) + (world_own - world_tgt)
    return shifted


def centroid_above(data: np.ndarray, threshold: float = 0.5) -> tuple[float, float] | None:
    """Sub-pixel centroid of all pixels above *threshold* (or None)."""
    from scipy.ndimage import center_of_mass

    mask = data > threshold
    if not mask.any():
        return None
    return center_of_mass(data, labels=mask, index=1)[::-1]


def make_exposures(*wcss: WCS):
    """Wrap WCS objects in objects exposing the ``.wcs``/``.shape`` contract."""
    from types import SimpleNamespace

    return [SimpleNamespace(wcs=w, shape=(64, 64)) for w in wcss]