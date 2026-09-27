"""Frame-center tables, footprint overlap statistics and overlap grouping."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy.wcs import WCS
from matplotlib.path import Path

from jwst_stack.io import CalExposure


@dataclass
class FrameInfo:
    """One row of the inspection table."""

    filename: str
    visit: int
    exposure: int
    ra_deg: float
    dec_deg: float
    pixscale_arcsec: float


def frame_table(exposures: list[CalExposure]) -> list[FrameInfo]:
    """Centers and pixel scale for every exposure, from the WCS."""
    rows: list[FrameInfo] = []
    for exp in exposures:
        ny, nx = exp.shape
        ra, dec = exp.wcs.all_pix2world(nx / 2.0, ny / 2.0, 0)
        scale = exp.wcs.proj_plane_pixel_scales()[0].to_value("arcsec")
        rows.append(
            FrameInfo(
                filename=exp.name,
                visit=exp.visit,
                exposure=exp.exposure,
                ra_deg=float(ra),
                dec_deg=float(dec),
                pixscale_arcsec=float(scale),
            )
        )
    return rows


def footprint_polygon(wcs: WCS, axes: tuple[int, int] = (2048, 2048)) -> np.ndarray:
    """World-coordinate footprint polygon corners (RA, Dec) of the frame.

    ``axes`` is the (naxis1, naxis2) pixel extent, i.e. ``shape[::-1]`` for a
    (ny, nx) array.
    """
    return wcs.calc_footprint(axes=axes)


def overlap_fraction(
    wcs_a: WCS,
    wcs_b: WCS,
    sample: int = 48,
    axes_a: tuple[int, int] = (2048, 2048),
    axes_b: tuple[int, int] = (2048, 2048),
) -> float:
    """Fraction of frame A's footprint overlapped by frame B.

    A coarse grid over frame A's sky bounding box is tested against both
    footprints via point-in-polygon.  ``sample`` is the number of grid cells
    along each axis; 48 gives milli-percent precision for grouping purposes.
    """
    pa = footprint_polygon(wcs_a, axes=axes_a)
    pb = footprint_polygon(wcs_b, axes=axes_b)
    lo_ra, hi_ra = pa[:, 0].min(), pa[:, 0].max()
    lo_dec, hi_dec = pa[:, 1].min(), pa[:, 1].max()
    ra = np.linspace(lo_ra, hi_ra, sample)
    dec = np.linspace(lo_dec, hi_dec, sample)
    grid = np.column_stack(
        (np.repeat(ra, sample), np.tile(dec, sample))
    )
    in_a = Path(pa).contains_points(grid)
    total = in_a.sum()
    if total == 0:
        return 0.0
    return float(np.logical_and(in_a, Path(pb).contains_points(grid)).sum() / total)


def pairwise_overlap(
    exposures: list[CalExposure], sample: int = 48
) -> np.ndarray:
    """Symmetric pair-wise overlap-fraction matrix for the exposures."""
    n = len(exposures)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            frac = overlap_fraction(
                exposures[i].wcs,
                exposures[j].wcs,
                sample=sample,
                axes_a=exposures[i].shape[::-1],
                axes_b=exposures[j].shape[::-1],
            )
            matrix[i, j] = frac
            matrix[j, i] = frac
    return matrix


def connected_groups(
    exposures: list[CalExposure],
    threshold: float = 0.15,
    sample: int = 48,
) -> tuple[list[list[int]], np.ndarray]:
    """Group indices whose footprints overlap by at least *threshold*.

    Returns (groups, matrix) where each group is a list of indices connected
    through pairwise overlap >= threshold.
    """
    matrix = pairwise_overlap(exposures, sample=sample)
    parent = list(range(len(exposures)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    n = len(exposures)
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i, j] >= threshold:
                union(i, j)

    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return [sorted(idx) for idx in groups.values()], matrix


def format_frame_table(rows: list[FrameInfo]) -> str:
    """Render the inspection table for the terminal."""
    head = f"{'file':<43} {'visit':>5} {'exp':>3} {'RA (deg)':>12} {'Dec (deg)':>12} {'arcsec/px':>9}"
    lines = [head]
    for r in sorted(rows, key=lambda x: (x.visit, x.exposure)):
        lines.append(
            f"{r.filename:<43} {r.visit:>5} {r.exposure:>3} "
            f"{r.ra_deg:12.7f} {r.dec_deg:12.7f} {r.pixscale_arcsec:9.5f}"
        )
    return "\n".join(lines)


def format_groups(
    exposures: list[CalExposure],
    groups: list[list[int]],
    matrix: np.ndarray,
    threshold: float,
) -> str:
    """Human-readable overlap grouping report."""
    lines = [f"Overlap threshold: >= {threshold:.0%} of a frame footprint"]
    for gi, group in enumerate(groups, start=1):
        member = ", ".join(
            f"v{e.visit}e{e.exposure}" for e in (exposures[i] for i in group)
        )
        lines.append(f"  group {gi} ({len(group)} frames): {member}")
    lines.append("\nPair-wise overlap fractions (rows/cols = group order):")
    head = " " * 9 + "".join(f"{i:>8}" for i in range(len(exposures)))
    lines.append(head)
    for i in range(len(exposures)):
        lines.append(
            f"{'['+str(exposures[i].visit)+'.'+str(exposures[i].exposure)+']':9}"
            + "".join(f"{matrix[i, j]:8.2f}" for j in range(len(exposures)))
        )
    return "\n".join(lines)