"""A single, fixed common output grid for the NGC 3324 mosaic.

The grid is derived once from MAST observation metadata: the union of the
``s_region`` footprint polygons of every NIRCam observation in proposal 2731
(all six filters).  It is written to ``out/grid.fits`` (WCS header + grid
shape) with a companion inputs cache (``out/grid_inputs/ngc3324_obs.csv``)
and is never rebuilt once it exists; every later pipeline stage loads this
file.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from jwst_stack.units import arcsec_to_deg

from jwst_stack import download

DEFAULT_PIXEL_SCALE_ARCSEC = 0.031
DEFAULT_PAD_PX = 30
FILTER_ORDER = ["F090W", "F187N", "F200W", "F335M", "F444W", "F444W;F470N"]


@dataclass
class Footprint:
    ra_min: float
    ra_max: float
    dec_min: float
    dec_max: float
    rows_used: int
    rows_total: int


def parse_s_region(value) -> np.ndarray | None:
    """Return (N, 2) polygon vertices (RA, Dec) for a MAST ``s_region`` string.

    Returns None for empty/masked values or non-POLYGON shapes.
    """
    if value is None or value is np.ma.masked:
        return None
    text = str(value).strip()
    if not text or not text.upper().startswith("POLYGON"):
        return None
    tokens = text.split()
    coords = tokens[1:]
    if len(coords) < 6 or len(coords) % 2 != 0:
        return None
    try:
        numbers = np.array([float(t) for t in coords])
    except ValueError:
        return None
    return numbers.reshape(-1, 2)


def union_footprint(rows: Table) -> Footprint:
    """Union bounding box over every usable MAST footprint row."""
    ra_lo, ra_hi = [], []
    dec_lo, dec_hi = [], []
    used = 0
    for row in rows:
        vertices = parse_s_region(row["s_region"])
        if vertices is None or len(vertices) == 0:
            continue
        used += 1
        ra_lo.append(vertices[:, 0].min())
        ra_hi.append(vertices[:, 0].max())
        dec_lo.append(vertices[:, 1].min())
        dec_hi.append(vertices[:, 1].max())
    if not ra_lo:
        raise ValueError("no usable s_region footprints in the MAST rows")
    return Footprint(
        ra_min=float(min(ra_lo)),
        ra_max=float(max(ra_hi)),
        dec_min=float(min(dec_lo)),
        dec_max=float(max(dec_hi)),
        rows_used=used,
        rows_total=len(rows),
    )


def build_grid_wcs(
    footprint: Footprint,
    pixel_scale_arcsec: float = DEFAULT_PIXEL_SCALE_ARCSEC,
    pad_px: int = DEFAULT_PAD_PX,
) -> tuple[WCS, tuple[int, int]]:
    """North-up TAN WCS covering *footprint* plus padding, and its shape."""
    scale_deg = arcsec_to_deg(pixel_scale_arcsec)
    nax1 = int(np.ceil((footprint.ra_max - footprint.ra_min) / scale_deg)) + 2 * pad_px
    nax2 = int(np.ceil((footprint.dec_max - footprint.dec_min) / scale_deg)) + 2 * pad_px
    shape = (nax2, nax1)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = ["deg", "deg"]
    wcs.wcs.crval = [
        0.5 * (footprint.ra_min + footprint.ra_max),
        0.5 * (footprint.dec_min + footprint.dec_max),
    ]
    wcs.wcs.crpix = [(nax1 + 1) / 2.0, (nax2 + 1) / 2.0]
    wcs.wcs.cdelt = [-scale_deg, scale_deg]
    return wcs, shape


def save_grid(
    path: str | Path,
    wcs: WCS,
    shape: tuple[int, int],
    pixel_scale_arcsec: float,
    footprint: Footprint,
    rows: Table,
) -> Path:
    path = Path(path)
    header = wcs.to_header(relax=True)
    header["GRDNX"] = (shape[1], "common grid naxis1 (pixels)")
    header["GRDNY"] = (shape[0], "common grid naxis2 (pixels)")
    header["GRDSCL"] = (pixel_scale_arcsec, "common grid pixel scale (arcsec)")
    header["GRDFRA1"] = (footprint.ra_min, "union footprint RA min (deg)")
    header["GRDFRA2"] = (footprint.ra_max, "union footprint RA max (deg)")
    header["GRDFDE1"] = (footprint.dec_min, "union footprint Dec min (deg)")
    header["GRDFDE2"] = (footprint.dec_max, "union footprint Dec max (deg)")
    header["GRDFRUS"] = (footprint.rows_used, "MAST rows with usable footprints")
    header["GRDFRTT"] = (footprint.rows_total, "MAST rows queried")
    header["HISTORY"] = "Common mosaic grid from proposal 2731 NIRCam footprints"
    header["HISTORY"] = "Built once from MAST observation metadata; never rebuilt"
    fits.PrimaryHDU(header=header).writeto(path, overwrite=True)
    return path


def load_grid(path: str | Path) -> tuple[WCS, tuple[int, int]]:
    """Read back (WCS, shape) for a grid produced by :func:`save_grid`."""
    path = Path(path)
    with fits.open(path) as hdul:
        header = hdul[0].header
    wcs = WCS(header)
    nax1 = int(header["GRDNX"])
    nax2 = int(header["GRDNY"])
    return wcs, (nax2, nax1)


def write_inputs_cache(rows: Table, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keep = rows[["obs_id", "filters", "s_region"]].copy()
    keep.write(str(path), format="ascii.csv", overwrite=True)
    return path


def read_inputs_cache(path: str | Path) -> Table | None:
    path = Path(path)
    if not path.exists():
        return None
    return Table.read(str(path), format="ascii.csv")


def fetch_obs_rows(
    proposal_id: int, cache_path: str | Path, refresh: bool = False
) -> Table:
    """MAST observation rows, preferring the on-disk cache when present."""
    if not refresh:
        cached = read_inputs_cache(cache_path)
        if cached is not None:
            return cached
    obs = download.query_obs_table(proposal_id=proposal_id)
    nircam = obs[[download.is_nircam_obs(o) for o in obs["obs_id"]]]
    write_inputs_cache(nircam, cache_path)
    return nircam


def summarize_grid(
    wcs: WCS,
    shape: tuple[int, int],
    footprints: Footprint,
    rows: Table,
    pixel_scale_arcsec: float,
) -> str:
    ny, nx = shape
    extent_deg_x = arcsec_to_deg(nx * pixel_scale_arcsec)
    extent_deg_y = arcsec_to_deg(ny * pixel_scale_arcsec)
    lines = [
        "common output grid (fixed, from MAST observation footprints)",
        f"  pixel scale : {pixel_scale_arcsec:.4f} arcsec/px",
        f"  shape       : {ny} x {nx} pixels (y, x)",
        f"  size        : {extent_deg_x * 60:.2f} x {extent_deg_y * 60:.2f} arcmin",
        f"  footprint RA : {footprints.ra_min:.6f}..{footprints.ra_max:.6f} deg",
        f"  footprint Dec: {footprints.dec_min:.6f}..{footprints.dec_max:.6f} deg",
        f"  rows used    : {footprints.rows_used}/{footprints.rows_total} MAST observations",
    ]
    by_filter = {}
    for row in rows:
        key = str(row["filters"])
        used = parse_s_region(row["s_region"]) is not None
        by_filter.setdefault(key, [0, 0])
        by_filter[key][0] += 1
        if not used:
            by_filter[key][1] += 1
    lines.append("  per-filter (MAST observation rows):")
    for key in sorted(by_filter):
        total, missing = by_filter[key]
        note = "" if missing == 0 else f"  ({missing} without usable footprint)"
        lines.append(f"    {key:<12} {total} obs {note}")
    return "\n".join(lines)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="jwst_stack grid")
    parser.add_argument("--outdir", default="out")
    parser.add_argument("--scale", type=float, default=DEFAULT_PIXEL_SCALE_ARCSEC)
    parser.add_argument("--pad", type=int, default=DEFAULT_PAD_PX)
    parser.add_argument("--proposal-id", type=int, default=download.PROPOSAL_ID)
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--show", action="store_true", help="print summary of an existing grid")
    return parser.parse_args(argv)


def run_grid(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    grid_path = outdir / "grid.fits"
    inputs_dir = outdir / "grid_inputs"
    cache_path = inputs_dir / "ngc3324_obs.csv"

    if grid_path.exists():
        wcs, shape = load_grid(grid_path)
        rows = read_inputs_cache(cache_path)
        if rows is None:
            print(f"grid exists but inputs cache missing at {cache_path}")
            return 1
        footprint = union_footprint(rows)
        scale = float(fits.getheader(grid_path, 0)["GRDSCL"])
        print(summarize_grid(wcs, shape, footprint, rows, scale))
        return 0

    print(f"querying MAST for proposal {args.proposal_id} NIRCam observations ...")
    rows = fetch_obs_rows(args.proposal_id, cache_path, refresh=args.refresh_cache)
    if len(rows) == 0:
        print("no NIRCam observations returned by MAST")
        return 1
    footprint = union_footprint(rows)

    wcs, shape = build_grid_wcs(footprint, args.scale, args.pad)
    outdir.mkdir(parents=True, exist_ok=True)
    save_grid(grid_path, wcs, shape, args.scale, footprint, rows)
    print(summarize_grid(wcs, shape, footprint, rows, args.scale))
    print(f"\nwrote {grid_path}")
    print(f"wrote inputs cache {cache_path} (grid is now fixed)")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    return run_grid(args)


if __name__ == "__main__":
    sys.exit(main())