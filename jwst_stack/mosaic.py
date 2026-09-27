"""Full-mosaic construction on the fixed common grid.

The Stage 4 mosaic spans every F200W exposure of the program (8 NIRCam short
wave detectors x 4 visits x 5 dithers = 160 frames) on the grid fixed in
``out/grid.fits`` (15895 x 22130 = 3.52e8 pixels).  A single float32 overlay of
that grid is 1.41 GB, so the existing in-memory :func:`jwst_stack.stack.
stack_visit` (which keeps every overlay of a group in RAM at once) would need
about 225 GB for 160 frames.  This module therefore streams the mosaic tile by
tile:

* registration is solved once per (visit, detector) group on a small crop of
  the *mosaic* grid, so the solutions are already expressed in mosaic-grid
  pixels and need no rescaling;
* each output tile is filled by reprojecting only the frames whose footprint
  intersects it, applying that frame's registration shift, and combining with
  the same sigma-clipping rule used per visit;
* the finished arrays are accumulated in disk-backed ``np.memmap`` scratch and
  streamed into a FITS file, so resident memory stays a few hundred MB
  regardless of mosaic size.

Interpolation follows the validated per-visit recipe: ``reproject_interp``
(bilinear) for the geometric resampling, then a cubic ``ndi_shift`` for the
sub-pixel registration translation with the NaN mask filled before and
restored after.  ``reproject`` must stay at bilinear order here: for
``order >= 2`` reproject 0.21 runs a *global* ``spline_filter`` over the input
image, which would propagate the NaNs of the DQ-masked pixels across the whole
frame (the same failure mode that makes ``ndi_shift`` need NaN filling).
"""

from __future__ import annotations

import json
import os
import time
from collections import OrderedDict
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from reproject import reproject_interp
from scipy.ndimage import shift as ndi_shift
from scipy.spatial import cKDTree

if os.name == "nt":
    import ctypes
    from ctypes import wintypes
else:  # pragma: no cover
    ctypes = None
    wintypes = None

from jwst_stack.grid import load_grid
from jwst_stack.io import CalExposure, mask_scaled_sci
from jwst_stack.register import (
    FrameShift,
    RegistrationResult,
    detect_stars,
    format_registration_summary,
    register_overlays,
)

DEFAULT_TILE_PX = 1024
DEFAULT_HALO_PX = 24
DEFAULT_CACHE_FRAMES = 32
_MIN_VALID_FOR_CLIP = 3
_NEGLIGIBLE_SHIFT_PX = 0.01


# --------------------------------------------------------------------------
# resources
# --------------------------------------------------------------------------


class _ProcessMemoryCounters(ctypes.Structure):
    """Windows ``PROCESS_MEMORY_COUNTERS`` (72 bytes on x64)."""

    _fields_ = [
        ("cb", wintypes.DWORD),
        ("PageFaultCount", wintypes.DWORD),
        ("PeakWorkingSetSize", ctypes.c_size_t),
        ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t),
        ("PeakPagefileUsage", ctypes.c_size_t),
    ]


def _process_memory() -> tuple[float, float]:
    """Return (peak working set, current working set) in GB.

    ``GetProcessMemoryInfo`` needs explicit argument types: without them ctypes
    truncates the 64-bit pseudo-handle and the call silently fails, which would
    have made every reported memory figure ``nan`` exactly when it matters.
    """
    if os.name != "nt":
        try:
            import resource

            peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            with open("/proc/self/statm", encoding="ascii") as fh:
                pages = int(fh.read().split()[1])
            return peak / 1e6, pages * os.sysconf("SC_PAGE_SIZE") / 1e9
        except Exception:
            return float("nan"), float("nan")
    try:
        counters = _ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(_ProcessMemoryCounters)
        handle = wintypes.HANDLE(ctypes.windll.kernel32.GetCurrentProcess())
        func = ctypes.windll.psapi.GetProcessMemoryInfo
        func.argtypes = [
            wintypes.HANDLE,
            ctypes.POINTER(_ProcessMemoryCounters),
            wintypes.DWORD,
        ]
        func.restype = wintypes.BOOL
        if func(handle, ctypes.byref(counters), counters.cb):
            return (
                counters.PeakWorkingSetSize / 1e9,
                counters.WorkingSetSize / 1e9,
            )
    except Exception:
        pass
    return float("nan"), float("nan")


def peak_rss_gb() -> float:
    """Peak working-set size of this process in GB."""
    return _process_memory()[0]


def current_rss_gb() -> float:
    """Current resident set size of this process in GB."""
    return _process_memory()[1]


def free_disk_gb(path: str | Path) -> float:
    """Free space on the volume holding *path*, in GB."""
    target = Path(path)
    while not target.exists() and target != target.parent:
        target = target.parent
    usage = __import__("shutil").disk_usage(target)
    return usage.free / 1e9


# --------------------------------------------------------------------------
# grid geometry
# --------------------------------------------------------------------------


def tile_sub_wcs(
    grid_wcs: WCS, y0: int, y1: int, x0: int, x1: int
) -> WCS:
    """WCS of the sub-array ``grid[y0:y1, x0:x1]`` of the fixed grid.

    Only ``CRPIX`` moves, so the sub-array shares the grid's projection,
    orientation, pixel scale and reference pixel world coordinate exactly.
    """
    sub = grid_wcs.deepcopy()
    sub.wcs.crpix = np.asarray(sub.wcs.crpix, dtype=float) - np.array(
        [x0, y0], dtype=float
    )
    sub.array_shape = (int(y1 - y0), int(x1 - x0))
    return sub


def _grid_bbox(
    exp: CalExposure, grid_wcs: WCS, shape: tuple[int, int], pad: int = 8
) -> tuple[int, int, int, int] | None:
    """Integer ``(y0, y1, x0, x1)`` bounding box of *exp* on the grid.

    The frame edges are sampled (not just the four corners) because the JWST
    ``SCI`` WCS carries SIP distortion and the detector is rotated relative to
    the north-up grid.  Returns None when the frame misses the grid entirely.
    """
    ny, nx = exp.shape
    xs = np.linspace(0.0, nx - 1.0, 17)
    ys = np.linspace(0.0, ny - 1.0, 17)
    edge_x, edge_y = np.meshgrid(xs, ys)
    mask = np.zeros((17, 17), dtype=bool)
    mask[0, :] = mask[-1, :] = True
    mask[:, 0] = mask[:, -1] = True
    px = edge_x[mask]
    py = edge_y[mask]
    ra, dec = exp.wcs.all_pix2world(px, py, 0)
    ra = np.atleast_1d(np.asarray(ra, dtype=float))
    dec = np.atleast_1d(np.asarray(dec, dtype=float))
    good = np.isfinite(ra) & np.isfinite(dec)
    if int(good.sum()) < 3:
        return None
    tangent = np.asarray(grid_wcs.wcs.crval, dtype=float)
    cos_sep = np.sin(np.radians(dec[good])) * np.sin(np.radians(tangent[1])) + np.cos(
        np.radians(dec[good])
    ) * np.cos(np.radians(tangent[1])) * np.cos(
        np.radians(ra[good] - tangent[0])
    )
    visible = cos_sep > 0.0
    if int(visible.sum()) < 3:
        return None
    gx, gy = grid_wcs.all_world2pix(ra[good][visible], dec[good][visible], 0)
    gx = np.atleast_1d(np.asarray(gx, dtype=float))
    gy = np.atleast_1d(np.asarray(gy, dtype=float))
    finite = np.isfinite(gx) & np.isfinite(gy)
    if int(finite.sum()) < 3:
        return None
    gx = gx[finite]
    gy = gy[finite]
    x0 = int(np.floor(gx.min())) - pad
    x1 = int(np.ceil(gx.max())) + 1 + pad
    y0 = int(np.floor(gy.min())) - pad
    y1 = int(np.ceil(gy.max())) + 1 + pad
    x0 = max(x0, 0)
    y0 = max(y0, 0)
    x1 = min(x1, int(shape[1]))
    y1 = min(y1, int(shape[0]))
    if x1 <= x0 or y1 <= y0:
        return None
    return y0, y1, x0, x1


def iter_tiles(
    shape: tuple[int, int], tile_px: int = DEFAULT_TILE_PX
) -> list[tuple[int, int, int, int]]:
    """Row-major ``(y0, y1, x0, x1)`` tiles covering the whole grid."""
    ny, nx = int(shape[0]), int(shape[1])
    tiles = []
    for y0 in range(0, ny, tile_px):
        for x0 in range(0, nx, tile_px):
            tiles.append((y0, min(y0 + tile_px, ny), x0, min(x0 + tile_px, nx)))
    return tiles


# --------------------------------------------------------------------------
# registration
# --------------------------------------------------------------------------


@dataclass
class FrameRegistration:
    """Registration solution for one frame, in mosaic-grid pixels."""

    filename: str
    visit: int
    exposure: int
    detector: str
    group: str
    dx: float
    dy: float
    background_offset: float
    n_detected: int
    n_matched: int
    residual_rms_px: float
    residual_max_px: float
    is_reference: bool = False
    visit_dx: float = 0.0
    visit_dy: float = 0.0

    @property
    def total_dx(self) -> float:
        """Intra-group shift plus the cross-visit correction, in grid pixels."""
        return self.dx + self.visit_dx

    @property
    def total_dy(self) -> float:
        return self.dy + self.visit_dy

    @property
    def applied(self) -> bool:
        return (
            max(
                abs(self.dx),
                abs(self.dy),
                abs(self.visit_dx),
                abs(self.visit_dy),
            )
            >= _NEGLIGIBLE_SHIFT_PX
        )


def group_by_visit_detector(
    exposures: list[CalExposure],
) -> list[tuple[tuple[int, str], list[int]]]:
    """Group frame indices by ``(visit, detector)``, ordered."""
    buckets: dict[tuple[int, str], list[int]] = {}
    for i, exp in enumerate(exposures):
        buckets.setdefault((exp.visit, exp.detector.lower()), []).append(i)
    return [
        (key, sorted(buckets[key], key=lambda i: (exposures[i].exposure, exposures[i].name)))
        for key in sorted(buckets)
    ]


def solve_registration(
    exposures: list[CalExposure],
    grid_wcs: WCS,
    grid_shape: tuple[int, int],
    fwhm_px: float = 3.0,
    match_radius_arcsec: float = 0.1,
    threshold_sigma: float = 5.0,
    bg_sigma: float = 3.0,
    interp_order: int = 3,
    cross_visit: bool = True,
    gauge_visit: int | None = None,
    verbose: bool = True,
) -> tuple[list[FrameRegistration], dict]:
    """Solve a translation + background offset per frame, group by group.

    Each ``(visit, detector)`` group is reprojected onto a small crop of the
    *fixed mosaic grid* and registered with the validated
    :func:`jwst_stack.register.register_overlays`.  Solving on a mosaic-grid
    crop (rather than on the per-visit grid) means the returned shifts need no
    pixel-scale conversion before they are applied in the mosaic pass.

    With ``cross_visit=True`` (the default) a second stage then solves one
    translation per visit from cross-visit star matches, because the per-group
    solve above cannot see inter-visit pointing error at all.
    """
    scale = float(grid_wcs.proj_plane_pixel_scales()[0].to_value("arcsec"))
    groups = group_by_visit_detector(exposures)
    solutions: list[FrameRegistration] = []
    reports: dict[str, dict] = {}

    for gi, ((visit, detector), members) in enumerate(groups, start=1):
        box = None
        for i in members:
            b = _grid_bbox(exposures[i], grid_wcs, grid_shape, pad=24)
            if b is None:
                continue
            if box is None:
                box = list(b)
            else:
                box[0] = min(box[0], b[0])
                box[1] = max(box[1], b[1])
                box[2] = min(box[2], b[2])
                box[3] = max(box[3], b[3])
        if box is None:
            for i in members:
                exp = exposures[i]
                solutions.append(
                    FrameRegistration(
                        filename=exp.name,
                        visit=exp.visit,
                        exposure=exp.exposure,
                        detector=exp.detector,
                        group=f"v{visit}_{detector}",
                        dx=0.0,
                        dy=0.0,
                        background_offset=0.0,
                        n_detected=0,
                        n_matched=0,
                        residual_rms_px=float("nan"),
                        residual_max_px=float("nan"),
                        is_reference=(i == members[0]),
                    )
                )
            continue
        y0, y1, x0, x1 = box
        sub_wcs = tile_sub_wcs(grid_wcs, y0, y1, x0, x1)
        overlays, footprints = [], []
        for i in members:
            exp = exposures[i]
            exp.load()
            cleaned = mask_scaled_sci(exp.sci, exp.dq)
            overlay, footprint = _reproject(cleaned, exp.wcs, sub_wcs)
            overlays.append(overlay)
            footprints.append(footprint)
            exp.release()
            del cleaned, overlay

        shifted, result = register_overlays(
            overlays,
            footprints,
            [exposures[i] for i in members],
            sub_wcs,
            pixel_scale_arcsec=scale,
            fwhm_px=fwhm_px,
            match_radius_arcsec=match_radius_arcsec,
            threshold_sigma=threshold_sigma,
            bg_sigma=bg_sigma,
            interp_order=interp_order,
        )
        reference_name = result.reference
        by_name = {f.filename: f for f in result.frames}
        for pos, i in enumerate(members):
            exp = exposures[i]
            if exp.name == reference_name:
                solutions.append(
                    FrameRegistration(
                        filename=exp.name,
                        visit=exp.visit,
                        exposure=exp.exposure,
                        detector=exp.detector,
                        group=f"v{visit}_{detector}",
                        dx=0.0,
                        dy=0.0,
                        background_offset=0.0,
                        n_detected=0,
                        n_matched=0,
                        residual_rms_px=float("nan"),
                        residual_max_px=float("nan"),
                        is_reference=True,
                    )
                )
                continue
            f = by_name[exp.name]
            solutions.append(
                FrameRegistration(
                    filename=exp.name,
                    visit=exp.visit,
                    exposure=exp.exposure,
                    detector=exp.detector,
                    group=f"v{visit}_{detector}",
                    dx=f.dx,
                    dy=f.dy,
                    background_offset=f.background_offset,
                    n_detected=f.n_detected,
                    n_matched=f.n_matched,
                    residual_rms_px=f.residual_rms_px,
                    residual_max_px=f.residual_max_px,
                    is_reference=False,
                )
            )
        reports[f"v{visit}_{detector}"] = {
            "reference": reference_name,
            "n_frames": len(members),
            "bbox": [y0, y1, x0, x1],
            "median_shift_px": result.median_shift_px,
            "mean_residual_rms_px": result.mean_residual_rms_px,
        }
        if verbose:
            shifts = [s for s in solutions if s.group == f"v{visit}_{detector}"]
            moved = [s for s in shifts if not s.is_reference]
            rms = [s.residual_rms_px for s in moved if np.isfinite(s.residual_rms_px)]
            med = (
                float(np.median([np.hypot(s.dx, s.dy) for s in moved])) if moved else 0.0
            )
            print(
                f"  [{gi:>2}/{len(groups)}] v{visit} {detector:<5} "
                f"{len(members)} frames  crop {y1-y0}x{x1-x0}  "
                f"median |shift| {med:.4f} px  mean residual rms "
                f"{np.mean(rms) if rms else float('nan'):.3f} px"
            )
        del overlays, footprints, shifted

    solutions.sort(key=lambda s: (s.visit, s.detector, s.exposure, s.filename))
    summary = {
        "pixel_scale_arcsec": scale,
        "n_frames": len(solutions),
        "n_groups": len(groups),
        "groups": reports,
    }
    if cross_visit:
        corrections, cross_report = solve_cross_visit_alignment(
            exposures,
            solutions,
            grid_wcs,
            fwhm_px=fwhm_px,
            threshold_sigma=threshold_sigma,
            gauge_visit=gauge_visit,
            verbose=verbose,
        )
        apply_cross_visit_corrections(solutions, corrections)
        summary["cross_visit"] = cross_report
    return solutions, summary


# --------------------------------------------------------------------------
# cross-visit alignment
# --------------------------------------------------------------------------


def _grid_star_positions(
    exposure: CalExposure,
    grid_wcs: WCS,
    fwhm_px: float,
    threshold_sigma: float,
    max_stars: int,
) -> np.ndarray:
    """Star positions of one frame expressed as mosaic-grid ``(x, y)`` pixels.

    The positions are routed through the frame WCS rather than through a
    reprojection so that the quantity measured is exactly the placement error
    that ``reproject_interp`` will later apply: the same SIP-distorted
    ``all_pix2world`` / ``all_world2pix`` pair.
    """
    exposure.load()
    cleaned = mask_scaled_sci(exposure.sci, exposure.dq)
    xs, ys, _flux = detect_stars(
        cleaned,
        fwhm_px=fwhm_px,
        threshold_sigma=threshold_sigma,
        max_stars=max_stars,
    )
    exposure.release()
    if xs.size == 0:
        return np.empty((0, 2), dtype=float)
    sky = exposure.wcs.all_pix2world(xs, ys, 0)
    gx, gy = grid_wcs.all_world2pix(sky[0], sky[1], 0)
    out = np.column_stack([np.asarray(gx, float), np.asarray(gy, float)])
    return out[np.isfinite(out).all(axis=1)]


def _match_catalogs(
    xy_a: np.ndarray, xy_b: np.ndarray, radius_px: float
) -> tuple[np.ndarray, np.ndarray] | None:
    """Mutual-nearest-neighbour match of two sky-position catalogues."""
    if xy_a.shape[0] < 1 or xy_b.shape[0] < 1:
        return None
    tree_b = cKDTree(xy_b)
    tree_a = cKDTree(xy_a)
    d_ab, j_ab = tree_b.query(xy_a)
    d_ba, i_ba = tree_a.query(xy_b)
    j_ab = np.asarray(j_ab, dtype=int)
    d_ab = np.asarray(d_ab, float)
    if j_ab.max(initial=0) >= xy_b.shape[0]:
        return None
    keep = d_ab <= radius_px
    if not np.any(keep):
        return None
    ia = np.flatnonzero(keep)
    ib = j_ab[ia]
    reciprocal = np.asarray(i_ba, dtype=int)[ib] == ia
    ia, ib = ia[reciprocal], ib[reciprocal]
    if ia.size == 0:
        return None
    return ia, ib


def _mad_sigma_clip(values: np.ndarray, clip_sigma: float = 3.0) -> np.ndarray:
    """Boolean mask of values within *clip_sigma* robust sigmas of the median."""
    if values.size == 0:
        return np.zeros(0, dtype=bool)
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma <= 0.0:
        return np.ones(values.size, dtype=bool)
    return np.abs(values - med) <= clip_sigma * sigma


def _solve_visit_translations(
    edges: list[dict], visits: list[int], gauge_visit: int | None = None
) -> tuple[dict[int, float], dict[int, float], dict]:
    """Weighted least squares for one translation per visit.

    Each edge constrains ``t[b] - t[a]`` to the measured offset of group *b*
    relative to group *a*.  The system is rank deficient by one, because a rigid
    translation of every frame is unobservable, so the solution must be gauged.

    The gauge is **not** the consensus of mutually-agreeing visits.  That was
    tried and is wrong: a visit can agree with the others yet be displaced in
    absolute terms, and gauging to the consensus silently redefines the
    mosaic's absolute sky frame.  Measured against the official i2d, that
    consensus frame sits 0.86 px from the reference, which is exactly the size
    of the inter-visit error being corrected - the correction was right and the
    origin was wrong.

    Instead the gauge is pinned to *gauge_visit*, which leaves that visit
    exactly where the header WCS put it and moves the others onto it.

    ``gauge_visit=None`` falls back to the lowest-numbered visit.  **That
    default is a heuristic, not a result, and it is recorded as such in
    ``gauge_source``.**  It is correct for F200W only because that filter's
    official i2d was measured to agree with visit 1's header WCS to 0.010 px.
    Which visit the official product is anchored to is a per-filter fact that
    must be re-measured, because a different filter can have a different visit
    displaced.  Pass ``gauge_visit`` explicitly once that measurement exists;
    internal self-consistency is not a substitute for it, since a mosaic can be
    perfectly self-consistent and still sit a full inter-visit error away from
    the reference.

    The consensus is still computed, but only as a diagnostic: it is what
    identifies a visit as displaced rather than merely mis-referenced.
    """
    index = {v: i for i, v in enumerate(visits)}
    n = len(edges)
    A = np.zeros((n, len(visits)))
    zx = np.zeros(n)
    zy = np.zeros(n)
    w = np.zeros(n)
    for e, edge in enumerate(edges):
        ia, ib = index[edge["visit_a"]], index[edge["visit_b"]]
        A[e, ia] = -1.0
        A[e, ib] = 1.0
        zx[e] = edge["dx"]
        zy[e] = edge["dy"]
        w[e] = np.sqrt(max(edge["n_matched"], 1)) / max(edge["mad_px"], 0.05)

    root_w = np.sqrt(w)
    Aw = A * root_w[:, None]
    tx, *_ = np.linalg.lstsq(Aw, zx * root_w, rcond=None)
    ty, *_ = np.linalg.lstsq(Aw, zy * root_w, rcond=None)

    res_x = A @ tx - zx
    res_y = A @ ty - zy
    resid_rms = float(np.sqrt(np.average(res_x**2 + res_y**2, weights=w) / 2.0))

    consensus = set(visits)
    floor_px = max(3.0 * resid_rms, 0.2)
    while len(consensus) > 3:
        means_x = float(np.mean([tx[index[v]] for v in consensus]))
        means_y = float(np.mean([ty[index[v]] for v in consensus]))
        worst = max(
            consensus,
            key=lambda v: np.hypot(tx[index[v]] - means_x, ty[index[v]] - means_y),
        )
        if np.hypot(tx[index[worst]] - means_x, ty[index[worst]] - means_y) <= floor_px:
            break
        consensus.discard(worst)

    if gauge_visit is None:
        gauge_visit = min(visits)
        gauge_source = "heuristic:min-visit"
    else:
        if gauge_visit not in index:
            raise ValueError(
                f"gauge_visit {gauge_visit} is not among the solved visits "
                f"{sorted(visits)}"
            )
        gauge_source = "explicit"
    tx = tx - tx[index[gauge_visit]]
    ty = ty - ty[index[gauge_visit]]

    return (
        {v: float(tx[index[v]]) for v in visits},
        {v: float(ty[index[v]]) for v in visits},
        {
            "residual_rms_px": resid_rms,
            "gauge_visit": gauge_visit,
            "gauge_source": gauge_source,
            "consensus_visits": sorted(consensus),
            "displaced_visits": sorted(set(visits) - consensus),
        },
    )


def solve_cross_visit_alignment(
    exposures: list[CalExposure],
    solutions: list[FrameRegistration],
    grid_wcs: WCS,
    fwhm_px: float = 3.0,
    threshold_sigma: float = 5.0,
    max_stars: int = 6000,
    match_radius_px: float = 10.0,
    min_matches: int = 20,
    gauge_visit: int | None = None,
    verbose: bool = True,
) -> tuple[dict[int, tuple[float, float]], dict]:
    """Solve one translation per *visit* from cross-visit star matches.

    ``solve_registration`` only ever compares frames inside a single
    ``(visit, detector)`` group, so inter-visit pointing error is invisible to
    it and passes straight into the mosaic.  This stage closes that gap: the
    reference frame of every group is star-detected in mosaic-grid coordinates,
    group pairs are cross-matched, and the resulting graph is solved for a
    per-visit translation.

    The solve is per *visit*, not per group, because the NIRCam short-wave
    detectors tile disjoint sky within a visit: the overlap graph contains no
    same-visit edges at all, so a per-detector term is unconstrained.  That in
    turn assumes pointing is constant across the 8 detectors of a visit, an
    assumption this dataset cannot test.
    """
    by_name = {e.name: e for e in exposures}
    references: list[tuple[tuple[int, str], str]] = []
    for rec in solutions:
        if rec.is_reference:
            references.append(((rec.visit, rec.detector), rec.filename))
    references.sort()

    catalogs: dict[tuple[int, str], np.ndarray] = {}
    for key, filename in references:
        exposure = by_name.get(filename)
        if exposure is None:
            continue
        catalogs[key] = _grid_star_positions(
            exposure, grid_wcs, fwhm_px, threshold_sigma, max_stars
        )
        if verbose:
            print(
                f"  cross-visit catalog v{key[0]} {key[1]:<5} "
                f"{catalogs[key].shape[0]:>5} stars"
            )

    keys = sorted(catalogs)
    edges: list[dict] = []
    for i, key_a in enumerate(keys):
        for key_b in keys[i + 1 :]:
            if key_a[0] == key_b[0]:
                continue
            matched = _match_catalogs(
                catalogs[key_a], catalogs[key_b], match_radius_px
            )
            if matched is None:
                continue
            ia, ib = matched
            if ia.size < min_matches:
                continue
            dx = catalogs[key_b][ib, 0] - catalogs[key_a][ia, 0]
            dy = catalogs[key_b][ib, 1] - catalogs[key_a][ia, 1]
            both = _mad_sigma_clip(dx) & _mad_sigma_clip(dy)
            if both.sum() < min_matches:
                continue
            dx, dy = dx[both], dy[both]
            med_x, med_y = float(np.median(dx)), float(np.median(dy))
            mad = float(
                np.median(np.abs(dx - med_x)) + np.median(np.abs(dy - med_y))
            )
            edges.append(
                {
                    "visit_a": key_a[0],
                    "detector_a": key_a[1],
                    "visit_b": key_b[0],
                    "detector_b": key_b[1],
                    "n_matched": int(both.sum()),
                    "dx": med_x,
                    "dy": med_y,
                    "mad_px": mad,
                }
            )

    visits = sorted({k[0] for k in keys})
    report: dict = {
        "n_groups": len(keys),
        "n_pairs": len(keys) * (len(keys) - 1) // 2,
        "n_edges": len(edges),
        "match_radius_px": match_radius_px,
        "min_matches": min_matches,
        "edges": edges,
        "per_detector_constancy_assumed": True,
    }
    if not edges:
        report["enabled"] = False
        report["reason"] = "no group pair reached the minimum matched-star count"
        if verbose:
            print("  cross-visit: no measurable edges; leaving header WCS")
        return {}, report

    tx, ty, fit = _solve_visit_translations(
        edges, visits, gauge_visit=gauge_visit
    )
    corrections = {v: (tx[v], ty[v]) for v in visits}
    report["enabled"] = True
    report["visits"] = {str(v): {"dx": tx[v], "dy": ty[v]} for v in visits}
    report.update(fit)
    if verbose:
        print(
            f"  cross-visit: {len(edges)} edges over {len(visits)} visits, "
            f"residual rms {fit['residual_rms_px']:.3f} px, "
            f"gauge v{fit['gauge_visit']} ({fit['gauge_source']})"
        )
        if fit["gauge_source"] != "explicit" and fit["displaced_visits"]:
            print(
                "    WARNING: gauge not verified against an external reference. "
                "Internal self-consistency does not pin the absolute frame; "
                "re-measure the raw frame WCSs against this filter's official "
                "i2d and pass --gauge-visit if the anchor is not the "
                "lowest-numbered visit."
            )
        for v in visits:
            print(f"    visit {v}: dx {tx[v]:+.4f}  dy {ty[v]:+.4f}  px")
    return corrections, report


def apply_cross_visit_corrections(
    solutions: list[FrameRegistration], corrections: dict[int, tuple[float, float]]
) -> None:
    """Fold per-visit corrections into every frame's solution, in place."""
    for rec in solutions:
        correction = corrections.get(rec.visit)
        if correction is None:
            continue
        rec.visit_dx, rec.visit_dy = float(correction[0]), float(correction[1])


def write_registration(
    solutions: list[FrameRegistration], summary: dict, path: str | Path
) -> Path:
    """Write the per-frame mosaic registration solution as JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(summary)
    payload["frames"] = [asdict(s) for s in solutions]
    moved = [s for s in solutions if not s.is_reference]
    payload["median_abs_shift_px"] = (
        float(np.median([np.hypot(s.dx, s.dy) for s in moved])) if moved else 0.0
    )
    rms = [s.residual_rms_px for s in moved if np.isfinite(s.residual_rms_px)]
    payload["mean_residual_rms_px"] = float(np.mean(rms)) if rms else 0.0
    payload["n_frames_shifted"] = int(sum(1 for s in moved if s.applied))
    payload["n_frames_wcs_kept"] = int(sum(1 for s in moved if not s.applied))
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def read_registration(path: str | Path) -> dict[str, FrameRegistration]:
    """Read a registration JSON written by :func:`write_registration`."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    out: dict[str, FrameRegistration] = {}
    for row in payload["frames"]:
        rec = FrameRegistration(**row)
        out[rec.filename] = rec
    return out


# --------------------------------------------------------------------------
# resampling and combination
# --------------------------------------------------------------------------


def _reproject(
    data: np.ndarray, input_wcs: WCS, output_wcs: WCS
) -> tuple[np.ndarray, np.ndarray]:
    """Bilinear reprojection onto *output_wcs*, NaN outside the footprint."""
    shape_out = tuple(int(s) for s in output_wcs.array_shape)
    result, footprint = reproject_interp(
        (data, input_wcs),
        output_projection=output_wcs,
        shape_out=shape_out,
        return_footprint=True,
    )
    result = np.array(result, copy=True)
    result[footprint <= 0] = np.nan
    return result, footprint


def apply_registration(
    overlay: np.ndarray, rec: FrameRegistration | None, interp_order: int = 3
) -> np.ndarray:
    """Apply a frame's registration solution to one reprojected overlay.

    Identical to the validated per-visit treatment: the NaN mask is filled
    before the cubic shift (a spline prefilter would otherwise spread the
    masked pixels across the whole tile) and restored afterwards, and the
    additive background offset is subtracted inside the footprint only.

    The intra-group shift and the cross-visit correction are summed into a
    single interpolation, so the NaN-mask discipline is exercised once.
    """
    if rec is None:
        return overlay
    valid = np.isfinite(overlay)
    bg = rec.background_offset
    if not rec.applied:
        if bg:
            return np.where(valid, overlay - bg, np.nan)
        return overlay
    moved = ndi_shift(
        np.where(valid, overlay, 0.0),
        shift=(-rec.total_dy, -rec.total_dx),
        order=interp_order,
        mode="nearest",
    )
    return np.where(valid, moved - bg, np.nan)


def combine_overlays(
    overlays: list[np.ndarray], sigma: float = 3.0, iterations: int = 3
) -> tuple[np.ndarray, np.ndarray]:
    """Sigma-clipped median combine of the tile overlays.

    Mirrors :func:`jwst_stack.stack.sigma_clip_stack`: pixels with at least
    three valid frames are clipped iteratively, thinner overlaps fall back to
    the plain median.  Returns (combined, coverage).
    """
    if not overlays:
        return np.empty(0), np.empty(0, dtype=np.int16)
    cube = np.stack(overlays, axis=0)
    valid = np.isfinite(cube)
    coverage = valid.sum(axis=0).astype(np.int16)
    clipable = coverage >= _MIN_VALID_FOR_CLIP
    if np.any(clipable):
        masked = sigma_clip(
            cube,
            sigma=sigma,
            maxiters=iterations,
            axis=0,
            masked=True,
            copy=True,
            cenfunc="median",
            stdfunc="mad_std",
        )
        clipped = masked.filled(np.nan)
        with np.errstate(invalid="ignore"):
            combined = np.nanmedian(clipped, axis=0)
    else:
        combined = np.full(cube.shape[1:], np.nan, dtype=np.float32)
    keep = ~clipable
    if np.any(keep):
        with np.errstate(invalid="ignore"):
            fallback = np.nanmedian(cube, axis=0)
        combined[keep] = fallback[keep]
    return combined, coverage


# --------------------------------------------------------------------------
# frame cache
# --------------------------------------------------------------------------


class _FrameCache:
    """Bounded LRU cache of DQ-masked frames keyed by filename."""

    def __init__(self, limit: int = DEFAULT_CACHE_FRAMES):
        self.limit = max(1, int(limit))
        self._items: "OrderedDict[str, tuple[np.ndarray, WCS]]" = OrderedDict()
        self.hits = 0
        self.loads = 0

    def get(self, exp: CalExposure) -> tuple[np.ndarray, WCS]:
        key = exp.name
        if key in self._items:
            self._items.move_to_end(key)
            self.hits += 1
            return self._items[key]
        exp.load()
        data = mask_scaled_sci(exp.sci, exp.dq)
        exp.release()
        self.loads += 1
        self._items[key] = (data, exp.wcs)
        while len(self._items) > self.limit:
            self._items.popitem(last=False)
        return self._items[key]


# --------------------------------------------------------------------------
# mosaic build
# --------------------------------------------------------------------------


@dataclass
class MosaicResult:
    """Summary of a completed mosaic build."""

    path: Path
    shape: tuple[int, int]
    n_frames: int
    n_tiles: int
    max_depth: int
    covered_pixels: int
    seconds: float
    peak_rss_gb: float
    scratch_gb: float
    output_gb: float
    frames_per_tile_max: int = 0
    cache_hits: int = 0
    cache_loads: int = 0

    def text(self) -> str:
        lines = [
            "mosaic build summary",
            f"  output          : {self.path}",
            f"  grid shape      : {self.shape[0]} x {self.shape[1]}",
            f"  frames combined : {self.n_frames}",
            f"  tiles           : {self.n_tiles} (max {self.frames_per_tile_max} frames/tile)",
            f"  covered pixels  : {self.covered_pixels:,} "
            f"({100.0 * self.covered_pixels / (self.shape[0] * self.shape[1]):.1f}% of grid)",
            f"  max depth       : {self.max_depth}",
            f"  wall time       : {self.seconds / 60.0:.1f} min",
            f"  peak RSS        : {self.peak_rss_gb:.2f} GB",
            f"  scratch (peak)  : {self.scratch_gb:.2f} GB",
            f"  output size     : {self.output_gb:.2f} GB",
            f"  frame loads     : {self.cache_loads} (cache hits {self.cache_hits})",
        ]
        return "\n".join(lines)


@dataclass
class MosaicPlan:
    """Pre-flight resource envelope for a mosaic build."""

    grid_shape: tuple[int, int]
    n_frames: int
    n_frames_in_grid: int
    n_tiles: int
    frames_per_tile_max: int
    mean_frames_per_tile: float
    reproject_calls: int
    scratch_gb: float
    output_gb: float
    cache_gb: float
    tile_peak_gb: float
    est_peak_rss_gb: float
    free_disk_gb: float
    frame_bytes_mb: float

    def text(self) -> str:
        n_px = self.grid_shape[0] * self.grid_shape[1]
        return "\n".join(
            [
                "Stage 4 mosaic pre-flight",
                f"  grid shape         : {self.grid_shape[0]} x {self.grid_shape[1]} "
                f"= {n_px:,} px ({n_px * 4 / 1e9:.2f} GB per float32 image)",
                f"  frames             : {self.n_frames} total, "
                f"{self.n_frames_in_grid} intersect the grid",
                f"  tiles              : {self.n_tiles} at {DEFAULT_TILE_PX} px",
                f"  frames per tile    : max {self.frames_per_tile_max}, "
                f"mean {self.mean_frames_per_tile:.1f}",
                f"  reproject calls    : {self.reproject_calls:,} "
                "(one per frame per overlapping tile, x2 passes if re-running)",
                "",
                "disk",
                f"  scratch (memmap)   : {self.scratch_gb:.2f} GB "
                "(SCI float32 + COVERAGE uint16, deleted after the FITS write)",
                f"  output FITS        : {self.output_gb:.2f} GB",
                f"  peak disk total    : {self.scratch_gb + self.output_gb:.2f} GB",
                f"  free disk          : {self.free_disk_gb:.1f} GB",
                "",
                "memory",
                f"  frame cache        : {self.cache_gb:.2f} GB "
                f"({int(self.cache_gb * 1e9 / (self.frame_bytes_mb * 1e6))} frames "
                f"x {self.frame_bytes_mb:.0f} MB)",
                f"  one tile's overlays: {self.tile_peak_gb:.2f} GB "
                "(worst tile, halo included)",
                f"  estimated peak RSS : {self.est_peak_rss_gb:.2f} GB "
                "(cache + overlays + sigma-clip temporaries)",
            ]
        )


def plan_mosaic(
    exposures: list[CalExposure],
    grid_path: str | Path,
    tile_px: int = DEFAULT_TILE_PX,
    halo_px: int = DEFAULT_HALO_PX,
    cache_frames: int = DEFAULT_CACHE_FRAMES,
    scratch_dir: str | Path | None = None,
) -> MosaicPlan:
    """Compute the disk and memory envelope of a mosaic build.

    Counted exactly from the per-frame bounding boxes, so the "frames per tile"
    and "reproject calls" figures are real rather than guessed.
    """
    grid_wcs, grid_shape = load_grid(grid_path)
    ny, nx = int(grid_shape[0]), int(grid_shape[1])
    n_px = ny * nx
    boxes = [_grid_bbox(exp, grid_wcs, grid_shape) for exp in exposures]
    usable = [b for b in boxes if b is not None]
    tiles = iter_tiles(grid_shape, tile_px)

    counts = []
    calls = 0
    for y0, y1, x0, x1 in tiles:
        ey0 = max(0, y0 - halo_px)
        ex0 = max(0, x0 - halo_px)
        ey1 = min(ny, y1 + halo_px)
        ex1 = min(nx, x1 + halo_px)
        k = 0
        for by0, by1, bx0, bx1 in usable:
            if by0 < ey1 and by1 > ey0 and bx0 < ex1 and bx1 > ex0:
                k += 1
        counts.append(k)
        calls += k

    frame_bytes = max((np.prod(exp.shape) * 4) for exp in exposures)
    tile_h = min(tile_px, ny) + 2 * halo_px
    tile_w = min(tile_px, nx) + 2 * halo_px
    frame_bytes_mb = frame_bytes / 1e6
    cache_gb = cache_frames * frame_bytes / 1e9
    tile_overlay_gb = max(counts) * tile_h * tile_w * 4 / 1e9
    clip_temporaries_gb = tile_overlay_gb * 3.5
    return MosaicPlan(
        grid_shape=(ny, nx),
        n_frames=len(exposures),
        n_frames_in_grid=len(usable),
        n_tiles=len(tiles),
        frames_per_tile_max=max(counts) if counts else 0,
        mean_frames_per_tile=(float(np.mean(counts)) if counts else 0.0),
        reproject_calls=calls,
        scratch_gb=n_px * 6 / 1e9,
        output_gb=n_px * 6 / 1e9,
        cache_gb=cache_gb,
        tile_peak_gb=tile_overlay_gb,
        est_peak_rss_gb=cache_gb + clip_temporaries_gb + 1.5,
        free_disk_gb=free_disk_gb(scratch_dir or Path.cwd()),
        frame_bytes_mb=frame_bytes_mb,
    )


def build_mosaic(
    exposures: list[CalExposure],
    grid_path: str | Path,
    out_path: str | Path,
    registrations: dict[str, FrameRegistration] | None = None,
    sigma: float = 3.0,
    iterations: int = 3,
    interp_order: int = 3,
    tile_px: int = DEFAULT_TILE_PX,
    halo_px: int = DEFAULT_HALO_PX,
    cache_frames: int = DEFAULT_CACHE_FRAMES,
    scratch_dir: str | Path | None = None,
    progress_every: int = 20,
) -> MosaicResult:
    """Build the full mosaic on the fixed grid, streaming tile by tile.

    *registrations* maps filename to the solution from :func:`solve_registration`
    (or :func:`read_registration`); frames without an entry keep their WCS
    position.  The output is a FITS file with a float32 ``SCI`` image and a
    ``COVERAGE`` depth image on the fixed grid WCS.
    """
    started = time.time()
    grid_wcs, grid_shape = load_grid(grid_path)
    ny, nx = int(grid_shape[0]), int(grid_shape[1])
    n_px = ny * nx

    boxes: list[tuple[int, int, int, int] | None] = [
        _grid_bbox(exp, grid_wcs, grid_shape) for exp in exposures
    ]
    usable = [i for i, b in enumerate(boxes) if b is not None]
    dropped = [exposures[i].name for i, b in enumerate(boxes) if b is None]
    if not usable:
        raise ValueError("no exposure footprint intersects the fixed grid")

    scratch = Path(scratch_dir) if scratch_dir else Path(out_path).parent
    scratch.mkdir(parents=True, exist_ok=True)
    sci_path = scratch / "mosaic_sci.f32"
    cov_path = scratch / "mosaic_cov.u16"
    n_bytes_sci = n_px * 4
    n_bytes_cov = n_px * 2
    scratch_gb = (n_bytes_sci + n_bytes_cov) / 1e9

    sci_mm = np.memmap(sci_path, dtype=np.float32, mode="w+", shape=(ny, nx))
    cov_mm = np.memmap(cov_path, dtype=np.uint16, mode="w+", shape=(ny, nx))
    cov_mm[:] = 0
    sci_mm[:] = np.nan

    cache = _FrameCache(cache_frames)
    tiles = iter_tiles(grid_shape, tile_px)
    max_depth = 0
    covered = 0
    frames_per_tile_max = 0

    print(
        f"mosaic: {len(usable)}/{len(exposures)} frames on grid, "
        f"{len(tiles)} tiles of {tile_px}px (halo {halo_px}px)"
    )
    if dropped:
        print(f"  skipped (outside grid): {len(dropped)}")
    print(
        f"  scratch {scratch_gb:.2f} GB, grid {ny}x{nx}, "
        f"free disk {free_disk_gb(scratch):.1f} GB"
    )

    for ti, (y0, y1, x0, x1) in enumerate(tiles, start=1):
        ey0 = max(0, y0 - halo_px)
        ex0 = max(0, x0 - halo_px)
        ey1 = min(ny, y1 + halo_px)
        ex1 = min(nx, x1 + halo_px)
        exp_wcs = tile_sub_wcs(grid_wcs, ey0, ey1, ex0, ex1)

        selected: list[int] = []
        for i in usable:
            by0, by1, bx0, bx1 = boxes[i]
            if by0 < ey1 and by1 > ey0 and bx0 < ex1 and bx1 > ex0:
                selected.append(i)
        frames_per_tile_max = max(frames_per_tile_max, len(selected))

        overlays: list[np.ndarray] = []
        for i in selected:
            exp = exposures[i]
            data, wcs = cache.get(exp)
            overlay, _ = _reproject(data, wcs, exp_wcs)
            moved = apply_registration(overlay, (registrations or {}).get(exp.name), interp_order)
            overlays.append(np.ascontiguousarray(moved[y0 - ey0 : y1 - ey0, x0 - ex0 : x1 - ex0]))
            del overlay, moved
        del exp_wcs

        if overlays:
            combined, coverage = combine_overlays(overlays, sigma, iterations)
            sci_mm[y0:y1, x0:x1] = combined
            cov_mm[y0:y1, x0:x1] = coverage
            if coverage.size:
                max_depth = max(max_depth, int(coverage.max()))
                covered += int(np.count_nonzero(coverage))
            del combined, coverage
        del overlays

        if progress_every and (ti % progress_every == 0 or ti == len(tiles)):
            el = time.time() - started
            rate = ti / el if el else 0.0
            print(
                f"  tile {ti:>4}/{len(tiles)}  {el / 60.0:5.1f} min  "
                f"{rate:5.2f} tile/s  rss {current_rss_gb():.2f} GB  "
                f"eta {(len(tiles) - ti) / rate / 60.0 if rate else float('nan'):5.1f} min",
                flush=True,
            )

    sci_mm.flush()
    cov_mm.flush()
    peak = peak_rss_gb()

    header = grid_wcs.to_header(relax=True)
    header["BUNIT"] = ("MJy/sr", "input SCI unit, preserved through stacking")
    header["FILTER"] = "F200W"
    header["NCOMBINE"] = (len(usable), "exposures combined")
    header["SIGCLIP"] = (sigma, "sigma-clip threshold")
    header["CLIPITER"] = (iterations, "sigma-clip iterations")
    header["NFRAME"] = (len(tiles), "output tiles")
    header["MAXDPTH"] = (max_depth, "maximum overlap depth")
    header["HISTORY"] = (
        f"Stage 4 full F200W mosaic: {len(usable)} exposures, "
        f"{len(set((e.visit, e.detector.lower()) for e in exposures))} "
        "visit x detector groups, tiled on the fixed common grid"
    )
    header["HISTORY"] = (
        "registration solved per visit x detector group on mosaic-grid crops; "
        f"shifts applied with scipy cubic shift, interp_order={interp_order}"
    )
    for i in usable:
        header["HISTORY"] = f"  {exposures[i].name}"

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print("  writing FITS from memmap ...", flush=True)
    with fits.HDUList(
        [
            fits.PrimaryHDU(data=sci_mm, header=header),
            fits.ImageHDU(data=cov_mm, name="COVERAGE"),
        ]
    ) as hdul:
        hdul.writeto(out_path, overwrite=True)
    output_gb = out_path.stat().st_size / 1e9
    peak = max(peak, peak_rss_gb())
    del sci_mm, cov_mm
    for p in (sci_path, cov_path):
        try:
            p.unlink()
        except OSError:
            pass

    result = MosaicResult(
        path=out_path,
        shape=(ny, nx),
        n_frames=len(usable),
        n_tiles=len(tiles),
        max_depth=max_depth,
        covered_pixels=covered,
        seconds=time.time() - started,
        peak_rss_gb=peak,
        scratch_gb=scratch_gb,
        output_gb=output_gb,
        frames_per_tile_max=frames_per_tile_max,
        cache_hits=cache.hits,
        cache_loads=cache.loads,
    )
    return result
