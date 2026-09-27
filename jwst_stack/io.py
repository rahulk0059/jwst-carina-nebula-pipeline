"""Discovery and lazy loading of JWST NIRCam Level 2 ``_cal.fits`` files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

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
    pupil: str = "CLEAR"
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

    ``pupil`` comes from the primary header.  It is not redundant with
    ``filter_name``: MAST's ``F444W;F470N`` is one *filter* exposed through
    two pupils, and ``FILTER`` reads ``F444W`` for both, so selecting on
    ``filter_name`` alone silently mixes two different bandpasses into one
    stack.  Every other filter in this program is ``PUPIL=CLEAR``.
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
        pupil=str(primary.get("PUPIL", "CLEAR")).strip().upper() or "CLEAR",
    )


def iter_exposures(input_dir: str | Path) -> Iterator[CalExposure]:
    """Yield ``CalExposure`` objects for every file, one at a time."""
    for path in find_cal_files(input_dir):
        yield read_cal_exposure(path)


def sky_pixel_scale_arcsec(exposure: CalExposure) -> float:
    """Pixel scale of an exposure in arcseconds."""
    scale = exposure.wcs.proj_plane_pixel_scales()[0]
    return float(scale.to_value("arcsec"))


def select_exposures(
    exposures: Iterable[CalExposure],
    detector: Sequence[str] | None = None,
    filter_name: Sequence[str] | None = None,
    pupil: Sequence[str] | None = None,
) -> list[CalExposure]:
    """Filter exposures by detector, filter and pupil, case-insensitively.

    ``filter_name`` matches ``FILTER`` only, which is why a bare ``F444W``
    request is ambiguous: MAST exposes F444W through both ``CLEAR`` and
    ``F470N``, and both read ``FILTER=F444W``.  Callers that pass ``filter_name``
    without ``pupil`` on a mixed set should check
    :func:`mixed_pupils` first rather than silently stacking two bandpasses.
    """
    def _norm(values: Sequence[str] | None) -> set[str] | None:
        if not values:
            return None
        return {str(v).strip().upper() for v in values if str(v).strip()}

    want_det = _norm(detector)
    want_filter = _norm(filter_name)
    want_pupil = _norm(pupil)
    out = []
    for e in exposures:
        if want_det and e.detector.strip().upper() not in want_det:
            continue
        if want_filter and e.filter_name.strip().upper() not in want_filter:
            continue
        if want_pupil and e.pupil.strip().upper() not in want_pupil:
            continue
        out.append(e)
    return out


def mixed_pupils(exposures: Iterable[CalExposure]) -> list[str]:
    """Sorted distinct pupils present in *exposures*; more than one is a hazard.

    Stacking two pupils of the same ``FILTER`` averages two different
    bandpasses, which is exactly the failure the long-wave ``F444W`` set invites.
    """
    return sorted({e.pupil.strip().upper() for e in exposures if e.pupil})


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