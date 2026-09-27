from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from jwst_stack import mosaic as M
from jwst_stack import stack as stack_mod
from jwst_stack.grid import load_grid
from jwst_stack.io import read_cal_exposure
from jwst_stack.validation import _match_tile_stars, validate_mosaic

from .conftest import make_wcs


# ---------------------------------------------------------------- geometry


def test_tile_sub_wcs_preserves_world_coordinates():
    grid = make_wcs(crpix=(100.5, 80.5))
    sub = M.tile_sub_wcs(grid, 20, 60, 10, 50)
    for x, y in [(10, 20), (25, 35), (49, 59), (10, 59)]:
        ra_g, dec_g = grid.all_pix2world(x, y, 0)
        ra_s, dec_s = sub.all_pix2world(x - 10, y - 20, 0)
        assert float(np.atleast_1d(ra_g)[0]) == pytest.approx(
            float(np.atleast_1d(ra_s)[0])
        )
        assert float(np.atleast_1d(dec_g)[0]) == pytest.approx(
            float(np.atleast_1d(dec_s)[0])
        )


def test_iter_tiles_covers_grid_exactly_once():
    shape = (250, 170)
    tiles = M.iter_tiles(shape, 100)
    seen = np.zeros(shape, dtype=np.int32)
    for y0, y1, x0, x1 in tiles:
        seen[y0:y1, x0:x1] += 1
    assert np.all(seen == 1)
    assert len(tiles) == 3 * 2


def test_grid_bbox_contains_frame_and_rejects_disjoint():
    from types import SimpleNamespace

    grid = make_wcs(crpix=(150.5, 150.5))
    shape = (300, 300)
    inside = SimpleNamespace(wcs=make_wcs(crpix=(150.5, 150.5)), shape=(64, 64))
    box = M._grid_bbox(inside, grid, shape)
    assert box is not None
    y0, y1, x0, x1 = box
    for px, py in ((0, 0), (63, 63), (0, 63), (63, 0), (32, 32)):
        gx, gy = grid.all_world2pix(*inside.wcs.all_pix2world(px, py, 0), 0)
        assert x0 - 0.5 <= float(np.atleast_1d(gx)[0]) <= x1 + 0.5
        assert y0 - 0.5 <= float(np.atleast_1d(gy)[0]) <= y1 + 0.5
    far = SimpleNamespace(
        wcs=make_wcs(crval=(10.0, 10.0), crpix=(32.5, 32.5)), shape=(64, 64)
    )
    assert M._grid_bbox(far, grid, shape) is None


# ------------------------------------------------------------ registration


def _star_image(cx, cy, amp=100.0, sig=2.0, shape=(200, 200)):
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    return amp * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sig**2))


def _centroid(data, cx, cy, half=24):
    from scipy.ndimage import center_of_mass

    y0, y1 = int(cy) - half, int(cy) + half
    x0, x1 = int(cx) - half, int(cx) + half
    sub = data[y0:y1, x0:x1].astype(float)
    cy_, cx_ = center_of_mass(sub)
    return float(cx_ + x0), float(cy_ + y0)


def _rec(dx, dy, bg=0.0):
    return M.FrameRegistration(
        filename="t.fits",
        visit=1,
        exposure=1,
        detector="NRCA1",
        group="g",
        dx=dx,
        dy=dy,
        background_offset=bg,
        n_detected=10,
        n_matched=5,
        residual_rms_px=0.1,
        residual_max_px=0.2,
    )


@pytest.mark.parametrize("dx,dy", [(0.37, -0.22), (-0.5, 0.5), (0.25, 0.25), (0.0, 0.0)])
def test_apply_registration_moves_star_by_negative_offset(dx, dy):
    """A recorded frame-minus-reference offset must be corrected by its negative."""
    img = _star_image(100.0, 100.0)
    cx, cy = _centroid(img, 100.0, 100.0)
    out = M.apply_registration(img, _rec(dx, dy), interp_order=3)
    nx, ny = _centroid(out, 100.0, 100.0)
    assert nx - cx == pytest.approx(-dx, abs=0.01)
    assert ny - cy == pytest.approx(-dy, abs=0.01)


def test_apply_registration_does_not_poison_nan():
    """The cubic spline prefilter must not spread the masked footprint."""
    img = _star_image(100.0, 100.0, shape=(220, 220))
    img[:, 150:] = np.nan
    out = M.apply_registration(img, _rec(0.3, -0.2), interp_order=3)
    assert np.isfinite(out[50:140, 40:140]).all()
    assert np.isnan(out[:, 150:]).all()


def test_apply_registration_skips_negligible_shift_but_subtracts_background():
    img = _star_image(100.0, 100.0) + 5.0
    out = M.apply_registration(img, _rec(0.0, 0.0, bg=2.0), interp_order=3)
    assert out[20, 20] == pytest.approx(3.0, abs=1e-5)
    assert out.max() == pytest.approx(img.max() - 2.0, rel=1e-6)


def test_apply_registration_keeps_nan_where_uncovered():
    img = np.full((120, 120), np.nan)
    img[50:70, 50:70] = 10.0
    out = M.apply_registration(img, _rec(0.2, 0.2, bg=1.0), interp_order=3)
    assert np.isnan(out[0, 0])
    assert np.isnan(out[100, 100])


# --------------------------------------------------------------- combining


def test_combine_overlays_matches_stack_semantics():
    rng = np.random.default_rng(3)
    base = rng.normal(10.0, 1.0, size=(40, 40))
    overlays = [base + rng.normal(0, 0.05, size=(40, 40)) for _ in range(5)]
    overlays[2] = overlays[2].copy()
    overlays[2][10:14, 10:14] += 500.0
    mine, cov = M.combine_overlays(overlays, sigma=3.0, iterations=3)
    theirs, cov2 = stack_mod.sigma_clip_stack(overlays, 3.0, 3, "median")
    assert np.allclose(np.nan_to_num(mine, nan=-1), np.nan_to_num(theirs, nan=-1))
    assert np.array_equal(cov, cov2)


def test_combine_overlays_handles_thin_overlap_without_clipping():
    overlays = [np.full((10, 10), np.nan) for _ in range(3)]
    overlays[0][2:5, 2:5] = 1.0
    overlays[1][2:5, 2:5] = 100.0
    combined, cov = M.combine_overlays(overlays)
    assert cov[3, 3] == 2
    assert np.isnan(combined[0, 0])
    assert np.isfinite(combined[3, 3])


# ------------------------------------------------------------- star match


def test_match_tile_stars_recovers_known_shift():
    rng = np.random.default_rng(11)
    ox = rng.uniform(0, 500, 60)
    oy = rng.uniform(0, 500, 60)
    dx, dy = 0.3, -0.2
    out = _match_tile_stars(ox, oy, ox + dx, oy + dy, radius_px=3.0)
    assert out["n_matched"] == 60
    assert out["dx_median"] == pytest.approx(dx, abs=0.05)
    assert out["dy_median"] == pytest.approx(dy, abs=0.05)
    assert out["nominal_frac"] == pytest.approx(1.0)


def test_match_tile_stars_rejects_far_matches():
    out = _match_tile_stars([0.0], [0.0], [50.0], [50.0], radius_px=2.0)
    assert out["n_matched"] == 0
    assert np.isnan(out["dx_median"])


# ------------------------------------------------------------ end-to-end


def _write_frame(path, wcs, data, visit=1, exposure=1, detector="NRCA1"):
    hdu = fits.PrimaryHDU()
    hdu.header["VISIT"] = visit
    hdu.header["EXPOSURE"] = exposure
    hdu.header["FILTER"] = "F200W"
    hdu.header["DETECTOR"] = detector
    hdu.header["EFFEXPTM"] = 100.0
    sci = fits.ImageHDU(data=np.asarray(data, dtype=np.float32), name="SCI")
    sci.header.update(wcs.to_header())
    err = fits.ImageHDU(data=np.ones_like(data, dtype=np.float32), name="ERR")
    dq = fits.ImageHDU(data=np.zeros(data.shape, dtype=np.uint32), name="DQ")
    fits.HDUList([hdu, sci, err, dq]).writeto(path, overwrite=True)


def _write_grid(path, wcs, shape):
    header = wcs.to_header(relax=True)
    header["GRDNX"] = (shape[1], "grid naxis1")
    header["GRDNY"] = (shape[0], "grid naxis2")
    fits.PrimaryHDU(header=header).writeto(path, overwrite=True)


def test_build_mosaic_end_to_end(tmp_path):
    """Two dithered frames of one group combine on the fixed grid."""
    frame_wcs = make_wcs(crpix=(32.5, 32.5))
    n = 64
    yy, xx = np.mgrid[0:n, 0:n]
    star = 100.0 * np.exp(-((xx - 32.0) ** 2 + (yy - 32.0) ** 2) / (2 * 2.0**2))
    shifted = 100.0 * np.exp(
        -((xx - 32.7) ** 2 + (yy - 31.6) ** 2) / (2 * 2.0**2)
    )
    f1 = tmp_path / "a.fits"
    f2 = tmp_path / "b.fits"
    _write_frame(f1, frame_wcs, star + 1.0, exposure=1)
    _write_frame(f2, frame_wcs, shifted + 1.0, exposure=2)
    exposures = [read_cal_exposure(f1), read_cal_exposure(f2)]

    grid_wcs = make_wcs(crpix=(48.5, 48.5))
    shape = (96, 96)
    grid_path = tmp_path / "grid.fits"
    _write_grid(grid_path, grid_wcs, shape)

    out = tmp_path / "mosaic.fits"
    result = M.build_mosaic(
        exposures, grid_path, out, registrations={}, tile_px=32, halo_px=8
    )
    assert result.shape == shape
    assert result.n_frames == 2
    assert out.exists()
    with fits.open(out) as hdul:
        sci = hdul[0].data
        cov = hdul["COVERAGE"].data
    assert sci.shape == shape
    assert int(cov.max()) == 2
    assert int(np.count_nonzero(cov)) > 0
    assert np.isfinite(sci[cov > 0]).all()
    assert result.output_gb > 0

    gx, gy = grid_wcs.all_world2pix(
        *frame_wcs.all_pix2world(32.0, 32.0, 0), 0
    )
    peak = np.unravel_index(np.nanargmax(sci), sci.shape)
    assert abs(peak[1] - float(np.atleast_1d(gx)[0])) <= 1.5
    assert abs(peak[0] - float(np.atleast_1d(gy)[0])) <= 1.5


def test_build_mosaic_applies_registration_solution(tmp_path):
    """A recorded offset must move the combined peak back onto the reference."""
    frame_wcs = make_wcs(crpix=(32.5, 32.5))
    n = 64
    yy, xx = np.mgrid[0:n, 0:n]
    star = 100.0 * np.exp(-((xx - 32.0) ** 2 + (yy - 32.0) ** 2) / (2 * 2.0**2))
    f1 = tmp_path / "a.fits"
    _write_frame(f1, frame_wcs, star + 1.0, exposure=1)
    exposures = [read_cal_exposure(f1)]
    grid_wcs = make_wcs(crpix=(48.5, 48.5))
    shape = (96, 96)
    grid_path = tmp_path / "grid.fits"
    _write_grid(grid_path, grid_wcs, shape)

    plain = tmp_path / "plain.fits"
    M.build_mosaic(exposures, grid_path, plain, registrations={}, tile_px=32, halo_px=8)
    with fits.open(plain) as hdul:
        base_peak = np.unravel_index(np.nanargmax(hdul[0].data), shape)

    shifted = tmp_path / "shifted.fits"
    M.build_mosaic(
        exposures,
        grid_path,
        shifted,
        registrations={"a.fits": _rec(0.5, -0.4)},
        tile_px=32,
        halo_px=8,
    )
    with fits.open(shifted) as hdul:
        new_peak = np.unravel_index(np.nanargmax(hdul[0].data), shape)
    assert new_peak[1] - base_peak[1] == pytest.approx(-0.5, abs=1)
    assert new_peak[0] - base_peak[0] == pytest.approx(0.4, abs=1)


def test_plan_mosaic_counts_frames_per_tile(tmp_path):
    from types import SimpleNamespace

    exposures = [
        SimpleNamespace(wcs=make_wcs(crpix=(100.5, 100.5)), shape=(64, 64)),
        SimpleNamespace(wcs=make_wcs(crpix=(140.5, 100.5)), shape=(64, 64)),
    ]
    grid_wcs = make_wcs(crpix=(128.5, 128.5))
    shape = (256, 256)
    grid_path = tmp_path / "grid.fits"
    _write_grid(grid_path, grid_wcs, shape)
    plan = M.plan_mosaic(exposures, grid_path, tile_px=64, halo_px=8)
    assert plan.grid_shape == shape
    assert plan.n_frames_in_grid == 2
    assert plan.frames_per_tile_max >= 1
    assert plan.reproject_calls >= 2
    assert plan.scratch_gb > 0
    assert "Stage 4 mosaic pre-flight" in plan.text()


def test_registration_json_roundtrip(tmp_path):
    rec = _rec(0.3, -0.2, bg=1.5)
    path = M.write_registration([rec], {"pixel_scale_arcsec": 0.031}, tmp_path / "r.json")
    back = M.read_registration(path)
    assert rec.filename in back
    assert back[rec.filename].dx == pytest.approx(0.3)
    assert back[rec.filename].dy == pytest.approx(-0.2)
    assert back[rec.filename].background_offset == pytest.approx(1.5)


def test_group_by_visit_detector_orders_groups():
    from types import SimpleNamespace

    def e(v, x, d, n):
        return SimpleNamespace(visit=v, exposure=x, detector=d, name=n)

    exposures = [
        e(2, 1, "NRCA2", "c"),
        e(1, 2, "NRCA1", "b"),
        e(1, 1, "NRCA1", "a"),
        e(1, 1, "nrcb1", "d"),
    ]
    groups = M.group_by_visit_detector(exposures)
    keys = [k for k, _ in groups]
    assert keys == [(1, "nrca1"), (1, "nrcb1"), (2, "nrca2")]
    assert [exposures[i].name for i in groups[0][1]] == ["a", "b"]


# ----------------------------------------------------------------- compare


def test_validate_mosaic_reports_offsets_and_noise(tmp_path):
    """A shifted, scaled copy of a mosaic must be measured as such."""
    grid_wcs = make_wcs(crpix=(48.5, 48.5))
    shape = (96, 96)
    grid_path = tmp_path / "grid.fits"
    _write_grid(grid_path, grid_wcs, shape)
    n = 96
    yy, xx = np.mgrid[0:n, 0:n]
    rng = np.random.default_rng(5)
    base = 10.0 + rng.normal(0, 0.1, size=(n, n))
    stars = np.zeros((n, n))
    for cx, cy in [(30, 30), (60, 40), (45, 70), (70, 20), (25, 60)]:
        stars += 100.0 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 1.5**2))
    truth = base + stars

    def write(path, data, wcs):
        hdu = fits.PrimaryHDU(data=np.asarray(data, dtype=np.float32))
        hdu.header.update(wcs.to_header())
        hdu.writeto(path, overwrite=True)

    mosaic_path = tmp_path / "mosaic.fits"
    i2d_path = tmp_path / "i2d.fits"
    write(mosaic_path, truth, grid_wcs)
    pix = 8.6e-6
    w2 = grid_wcs.deepcopy()
    w2.wcs.crval = np.asarray(w2.wcs.crval) + np.array([-0.5, 0.4]) * pix
    write(i2d_path, truth + 0.25, w2)

    res = validate_mosaic(
        mosaic_path,
        i2d_path,
        outdir=tmp_path,
        tile_px=48,
        halo_px=8,
        write_diff=False,
        progress_every=0,
    )
    assert res.overlap_pixels > 0
    assert res.background_difference == pytest.approx(-0.25, abs=0.15)
    assert res.clean_mad == pytest.approx(0.0, abs=0.2)
    assert res.matched_stars >= 3
    assert abs(res.star_dx_median) < 1.5
    assert abs(res.star_dy_median) < 1.5
    assert (tmp_path / f"{mosaic_path.stem}_vs_{i2d_path.stem}.json").exists()


# ---------------------------------------------------------- cross-visit


def _star_field(cols=10, rows=6, pitch=12.0, seed=0, origin=(20.0, 20.0)):
    """A jittered grid of star positions with a guaranteed large separation."""
    rng = np.random.default_rng(seed)
    pts = [
        (
            origin[0] + c * pitch + rng.uniform(-1.0, 1.0),
            origin[1] + r * pitch + rng.uniform(-1.0, 1.0),
        )
        for r in range(rows)
        for c in range(cols)
    ]
    return np.array(pts, dtype=float)


def _reference_solution(visit, detector):
    return M.FrameRegistration(
        filename=f"v{visit}_{detector}_ref.fits",
        visit=visit,
        exposure=1,
        detector=detector,
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


def _group_plan(visits=(1, 2), detectors=("nrca1", "nrcb1")):
    """Exposure stand-ins plus one reference registration per group."""
    from types import SimpleNamespace

    exposures, solutions = [], []
    for visit in visits:
        for detector in detectors:
            name = f"v{visit}_{detector}_ref.fits"
            exposures.append(
                SimpleNamespace(
                    name=name,
                    visit=visit,
                    exposure=1,
                    detector=detector,
                    wcs=make_wcs(crpix=(32.5, 32.5)),
                )
            )
            solutions.append(_reference_solution(visit, detector))
    return exposures, solutions


def test_match_catalogs_recovers_known_offset():
    field = _star_field()
    shifted = field + (0.6, -0.4)
    ia, ib = M._match_catalogs(field, shifted, 5.0)
    assert ia is not None
    assert ia.size == field.shape[0]
    dx = shifted[ib, 0] - field[ia, 0]
    dy = shifted[ib, 1] - field[ia, 1]
    assert np.median(dx) == pytest.approx(0.6, abs=0.01)
    assert np.median(dy) == pytest.approx(-0.4, abs=0.01)


def test_match_catalogs_rejects_offset_outside_radius():
    field = _star_field()
    far = field + (500.0, 0.0)
    assert M._match_catalogs(field, far, 5.0) is None


def test_solve_visit_translations_gauges_to_anchor_visit():
    """The measured inter-visit error: v1 displaced, v2/3/4 in agreement.

    The gauge is the anchor visit, *not* the consensus: visit 1 is the visit
    the official i2d agrees with, so pinning it to zero preserves the absolute
    sky frame while still removing the relative misalignment.
    """
    truth = {1: (-0.49, -0.65), 2: (0.12, 0.04), 3: (0.07, 0.23), 4: (0.18, 0.22)}
    edges = [
        {
            "visit_a": a,
            "detector_a": f"nrca{k}",
            "visit_b": b,
            "detector_b": f"nrcb{k}",
            "n_matched": 500,
            "dx": truth[b][0] - truth[a][0],
            "dy": truth[b][1] - truth[a][1],
            "mad_px": 0.06,
        }
        for a, b in ((2, 1), (3, 1), (4, 1))
        for k in (1, 2)
    ]
    tx, ty, fit = M._solve_visit_translations(edges, [1, 2, 3, 4])
    assert fit["gauge_visit"] == 1
    assert fit["displaced_visits"] == [1]
    assert set(fit["consensus_visits"]) == {2, 3, 4}
    # The anchor visit is left exactly where the header WCS put it.
    assert tx[1] == pytest.approx(0.0, abs=1e-12)
    assert ty[1] == pytest.approx(0.0, abs=1e-12)
    for v in (2, 3, 4):
        assert tx[v] == pytest.approx(truth[v][0] - truth[1][0], abs=0.01)
        assert ty[v] == pytest.approx(truth[v][1] - truth[1][1], abs=0.01)
    assert fit["residual_rms_px"] < 1e-6


def test_solve_visit_translations_honours_an_explicit_gauge():
    """Re-gauging must only shift the solution, never change the relative fit."""
    truth = {1: (-0.49, -0.65), 2: (0.12, 0.04), 3: (0.07, 0.23), 4: (0.18, 0.22)}
    edges = [
        {
            "visit_a": a,
            "detector_a": "d",
            "visit_b": b,
            "detector_b": "e",
            "n_matched": 200,
            "dx": truth[b][0] - truth[a][0],
            "dy": truth[b][1] - truth[a][1],
            "mad_px": 0.06,
        }
        for a, b in ((1, 2), (2, 3), (3, 4))
    ]
    tx1, ty1, _ = M._solve_visit_translations(edges, [1, 2, 3, 4], gauge_visit=1)
    tx3, ty3, fit3 = M._solve_visit_translations(edges, [1, 2, 3, 4], gauge_visit=3)
    assert fit3["gauge_visit"] == 3
    assert tx3[3] == pytest.approx(0.0, abs=1e-12)
    for v in (1, 2, 3, 4):
        assert tx3[v] - tx1[v] == pytest.approx(-(tx1[3] - 0.0), abs=1e-9)
        assert ty3[v] - ty1[v] == pytest.approx(-(ty1[3] - 0.0), abs=1e-9)
    # Relative offsets are gauge-invariant.
    for v in (1, 2, 4):
        assert tx3[v] - tx3[3] == pytest.approx(tx1[v] - tx1[3], abs=1e-9)


def _chain_edges(truth, pairs):
    return [
        {
            "visit_a": a,
            "detector_a": "d",
            "visit_b": b,
            "detector_b": "e",
            "n_matched": 200,
            "dx": truth[b][0] - truth[a][0],
            "dy": truth[b][1] - truth[a][1],
            "mad_px": 0.06,
        }
        for a, b in pairs
    ]


def test_gauge_source_records_whether_the_anchor_was_measured_or_guessed():
    """A guessed gauge must be distinguishable from a verified one.

    The min-visit default is a heuristic that happens to be right for F200W.
    A different filter can have a different visit displaced, and internal
    self-consistency cannot detect that, so the report must not present a guess
    and a measurement identically.
    """
    truth = {1: (-0.49, -0.65), 2: (0.12, 0.04), 3: (0.07, 0.23), 4: (0.18, 0.22)}
    edges = _chain_edges(truth, ((1, 2), (2, 3), (3, 4)))

    _, _, guessed = M._solve_visit_translations(edges, [1, 2, 3, 4])
    assert guessed["gauge_visit"] == 1
    assert guessed["gauge_source"] == "heuristic:min-visit"

    _, _, explicit = M._solve_visit_translations(
        edges, [1, 2, 3, 4], gauge_visit=2
    )
    assert explicit["gauge_visit"] == 2
    assert explicit["gauge_source"] == "explicit"


def test_solve_visit_translations_rejects_a_gauge_visit_outside_the_solve():
    truth = {1: (0.0, 0.0), 2: (0.1, 0.1)}
    edges = _chain_edges(truth, ((1, 2),))
    with pytest.raises(ValueError, match="not among the solved visits"):
        M._solve_visit_translations(edges, [1, 2], gauge_visit=7)


def test_ungauged_cross_visit_solve_warns_even_when_internally_consistent():
    """A displaced visit must not be reported as if it were anchored.

    This is the exact shape of the original bug: the solve is internally
    self-consistent to the noise floor, and only an external reference reveals
    that the frame is in the wrong place.  The heuristic gauge must therefore
    announce itself.
    """
    truth = {1: (0.9, 0.9), 2: (0.1, 0.1), 3: (0.12, 0.09), 4: (0.08, 0.11)}
    edges = _chain_edges(truth, ((1, 2), (2, 3), (3, 4)))

    # Visit 1 is the displaced one here, so the min-visit heuristic happens to
    # be right -- the point is that nothing internal can tell the difference,
    # and gauge_source must not pretend otherwise.
    _, _, fit = M._solve_visit_translations(edges, [1, 2, 3, 4])
    assert fit["displaced_visits"] == [1]
    assert fit["gauge_source"] == "heuristic:min-visit"

    # Now the same geometry with a *different* visit displaced: the solve is
    # just as self-consistent, the residual is still at the noise floor, and
    # only gauge_source distinguishes this from a verified anchor.
    truth2 = {1: (0.1, 0.1), 2: (0.11, 0.09), 3: (0.12, 0.1), 4: (0.95, 0.9)}
    edges2 = _chain_edges(truth2, ((1, 2), (2, 3), (3, 4)))
    _, _, fit2 = M._solve_visit_translations(edges2, [1, 2, 3, 4])
    assert fit2["displaced_visits"] == [4]
    assert fit2["residual_rms_px"] == pytest.approx(fit["residual_rms_px"], rel=1e-6)
    assert fit2["gauge_visit"] == 1
    assert fit2["gauge_source"] == "heuristic:min-visit"
    # Gauging to the real anchor recovers it, and is the only way to know.
    _, _, fit2b = M._solve_visit_translations(edges2, [1, 2, 3, 4], gauge_visit=4)
    assert fit2b["gauge_source"] == "explicit"
    assert fit2b["gauge_visit"] == 4


def test_cli_parses_gauge_visit():
    import jwst_stack.cli as cli

    args = cli.build_parser().parse_args(
        [
            "mosaic",
            "--detector", "nrca1",
            "--filter", "F200W",
            "--grid", "g.fits",
            "--out", "m.fits",
            "--registration", "reg.json",
            "--gauge-visit", "2",
        ]
    )
    assert args.gauge_visit == 2
    assert args.recompute_registration is False
    assert args.no_cross_visit is False
    default = cli.build_parser().parse_args(
        [
            "mosaic",
            "--detector", "nrca1",
            "--filter", "F200W",
            "--grid", "g.fits",
            "--out", "m.fits",
            "--registration", "reg.json",
        ]
    )
    assert default.gauge_visit is None


def test_solve_visit_translations_identifies_a_displaced_visit():
    truth = {1: (-0.10, 0.05), 2: (0.02, -0.01), 3: (0.03, 0.04), 4: (-0.01, 0.02)}
    edges = [
        {
            "visit_a": a,
            "detector_a": "d",
            "visit_b": b,
            "detector_b": "e",
            "n_matched": 100,
            "dx": truth[b][0] - truth[a][0],
            "dy": truth[b][1] - truth[a][1],
            "mad_px": 0.05,
        }
        for a, b in ((1, 2), (2, 3), (3, 4))
    ]
    tx, ty, fit = M._solve_visit_translations(edges, [1, 2, 3, 4])
    # All four agree, so none is displaced and the gauge is just the anchor.
    assert fit["displaced_visits"] == []
    assert set(fit["consensus_visits"]) == {1, 2, 3, 4}
    assert tx[1] == pytest.approx(0.0, abs=1e-12)
    for v in (2, 3, 4):
        assert tx[v] == pytest.approx(truth[v][0] - truth[1][0], abs=0.01)


def test_cross_visit_alignment_recovers_known_visit_offset(monkeypatch):
    """A pointing error confined to visit 1 must be recovered up to the gauge."""
    exposures, solutions = _group_plan()
    field = _star_field()
    truth = {1: (0.80, -0.60), 2: (0.0, 0.0)}
    monkeypatch.setattr(
        M,
        "_grid_star_positions",
        lambda exp, *args, **kwargs: field + truth[exp.visit],
    )
    corrections, report = M.solve_cross_visit_alignment(
        exposures, solutions, make_wcs(crpix=(100.5, 100.5)), verbose=False
    )
    assert report["enabled"] is True
    assert report["n_edges"] == 4
    assert report["residual_rms_px"] < 0.05
    # Only the difference between visits is observable; the gauge splits it.
    assert corrections[1][0] - corrections[2][0] == pytest.approx(0.80, abs=0.02)
    assert corrections[1][1] - corrections[2][1] == pytest.approx(-0.60, abs=0.02)


def test_cross_visit_alignment_is_disabled_without_edges(monkeypatch):
    exposures, solutions = _group_plan(visits=(1,))
    field = _star_field()
    monkeypatch.setattr(
        M, "_grid_star_positions", lambda exp, *args, **kwargs: field
    )
    corrections, report = M.solve_cross_visit_alignment(
        exposures, solutions, make_wcs(crpix=(100.5, 100.5)), verbose=False
    )
    assert corrections == {}
    assert report["enabled"] is False
    assert report["n_edges"] == 0


def test_cross_visit_alignment_flags_per_detector_constancy_assumption(monkeypatch):
    exposures, solutions = _group_plan()
    field = _star_field()
    monkeypatch.setattr(M, "_grid_star_positions", lambda exp, *a, **k: field)
    _corrections, report = M.solve_cross_visit_alignment(
        exposures, solutions, make_wcs(crpix=(100.5, 100.5)), verbose=False
    )
    assert report["per_detector_constancy_assumed"] is True


def test_apply_cross_visit_corrections_moves_even_a_reference_frame():
    """Reference frames carry no intra-group shift but must still be corrected."""
    rec = _rec(0.0, 0.0)
    rec.is_reference = True
    assert rec.applied is False
    M.apply_cross_visit_corrections([rec], {1: (0.5, -0.5)})
    assert rec.applied is True
    assert rec.visit_dx == pytest.approx(0.5)
    assert rec.visit_dy == pytest.approx(-0.5)
    assert rec.total_dx == pytest.approx(0.5)


def test_apply_registration_sums_intra_and_cross_visit_shifts():
    img = _star_image(100.0, 100.0)
    rec = _rec(0.2, -0.1)
    rec.visit_dx, rec.visit_dy = 0.3, 0.4
    assert rec.total_dx == pytest.approx(0.5)
    assert rec.total_dy == pytest.approx(0.3)
    cx, cy = _centroid(img, 100.0, 100.0)
    out = M.apply_registration(img, rec, interp_order=3)
    nx, ny = _centroid(out, 100.0, 100.0)
    assert nx - cx == pytest.approx(-0.5, abs=0.01)
    assert ny - cy == pytest.approx(-0.3, abs=0.01)


def test_apply_registration_sums_shifts_without_poisoning_nan():
    img = _star_image(100.0, 100.0, shape=(220, 220))
    img[:, 150:] = np.nan
    rec = _rec(0.3, -0.2)
    rec.visit_dx, rec.visit_dy = 0.2, 0.15
    out = M.apply_registration(img, rec, interp_order=3)
    assert np.isfinite(out[50:140, 40:140]).all()
    assert np.isnan(out[:, 150:]).all()


def test_registration_json_roundtrip_preserves_cross_visit(tmp_path):
    rec = _rec(0.3, -0.2, bg=1.5)
    rec.visit_dx, rec.visit_dy = -0.49, -0.65
    path = M.write_registration(
        [rec], {"pixel_scale_arcsec": 0.031}, tmp_path / "r.json"
    )
    back = M.read_registration(path)[rec.filename]
    assert back.visit_dx == pytest.approx(-0.49)
    assert back.visit_dy == pytest.approx(-0.65)
    assert back.total_dx == pytest.approx(0.3 - 0.49)
    assert back.total_dy == pytest.approx(-0.2 - 0.65)


def test_registration_json_without_visit_fields_defaults_to_zero(tmp_path):
    """Registration JSONs written before the cross-visit stage must still load."""
    import json

    rec = _rec(0.3, -0.2)
    path = M.write_registration(
        [rec], {"pixel_scale_arcsec": 0.031}, tmp_path / "r.json"
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    for row in payload["frames"]:
        row.pop("visit_dx")
        row.pop("visit_dy")
    path.write_text(json.dumps(payload), encoding="utf-8")
    back = M.read_registration(path)[rec.filename]
    assert back.visit_dx == 0.0
    assert back.visit_dy == 0.0
    assert back.dx == pytest.approx(0.3)
