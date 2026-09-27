"""Discovery and lazy loading of JWST NIRCam Level 2 ``_cal.fits`` files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

_DO_NOT_USE_BIT = 1 << 0


@dataclass
class CalExposure:
    """Metadata plus on-demand data access for a single calibrated exposure.

    Science arrays are kept on disk until :meth:`load` is called; only one
    exposure is held in memory at a time by the stacking pipeline.
    """

    path: Path
    visit: int
    exposure: int
    filter_name: str
    detector: str
    effexptm: float
    wcs: WCS
    shape: tuple[int, int]
    _sci: np.ndarray | None = None
    _err: np.ndarray | None = None
    _dq: np.ndarray | None = None

    @property
    def name(self) -> str:
        return self.path.name

    @property
    def sci(self) -> np.ndarray:
        if self._sci is None:
            self.load()
        return self._sci

    @property
    def err(self) -> np.ndarray:
        if self._err is None:
            self.load()
        return self._err

    @property
    def dq(self) -> np.ndarray:
        if self._dq is None:
            self.load()
        return self._dq

    def load(self) -> None:
        """Read SCI/ERR/DQ into memory.

        The DQ extension cannot be memory-mapped because it carries
        BZERO/BSCALE/BLANK scaling, so the file is opened without memmap.
        """
        with fits.open(self.path) as hdul:
            sci = hdul["SCI"].data
            err = hdul["ERR"].data
            dq = hdul["DQ"].data.astype(np.int64)
        self._sci = np.asarray(sci)
        self._err = np.asarray(err)
        self._dq = dq

    def release(self) -> None:
        """Drop loaded arrays to free memory."""
        self._sci = None
        self._err = None
        self._dq = None


def find_cal_files(input_dir: str | Path) -> list[Path]:
    """Return every ``*_cal.fits`` file below *input_dir*, sorted."""
    root = Path(input_dir)
    return sorted(p for p in root.rglob("*_cal.fits") if p.is_file())


def read_cal_exposure(path: str | Path) -> CalExposure:
    """Build a :class:`CalExposure` from a Level 2 file, loading only headers.

    The WCS is read from the SCI extension header (RA---TAN-SIP), which
    ``reproject`` consumes directly.
    """
    path = Path(path)
    with fits.open(path, memmap=False) as hdul:
        primary = hdul[0].header
        sci_header = hdul["SCI"].header
    wcs = WCS(sci_header)
    ny = sci_header["NAXIS2"]
    nx = sci_header["NAXIS1"]
    return CalExposure(
        path=path,
        visit=int(primary.get("VISIT")),
        exposure=int(primary.get("EXPOSURE")),
        filter_name=str(primary.get("FILTER", "")),
        detector=str(primary.get("DETECTOR", "")),
        effexptm=float(primary.get("EFFEXPTM", 1.0)),
        wcs=wcs,
        shape=(ny, nx),
    )


def iter_exposures(input_dir: str | Path) -> Iterator[CalExposure]:
    """Yield ``CalExposure`` objects for every file, one at a time."""
    for path in find_cal_files(input_dir):
        yield read_cal_exposure(path)


def sky_pixel_scale_arcsec(exposure: CalExposure) -> float:
    """Pixel scale of an exposure in arcseconds."""
    scale = exposure.wcs.proj_plane_pixel_scales()[0]
    return float(scale.to_value("arcsec"))


def mask_scaled_sci(sci: np.ndarray, dq: np.ndarray | None) -> np.ndarray:
    """Return a NaN-masked copy of the SCI array.

    Pixels are set to NaN when the DQ ``DO_NOT_USE`` bit is set, or when the
    value is not finite.  One percent of real NIRCam pixels are flagged and
    already NaN; this keeps all of them out of the stack.
    """
    out = np.array(sci, dtype=np.float32, copy=True)
    bad = ~np.isfinite(out)
    if dq is not None:
        bad |= (dq & _DO_NOT_USE_BIT) != 0
    out[bad] = np.nan
    return out