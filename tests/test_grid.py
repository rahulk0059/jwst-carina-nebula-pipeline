"""Tests for the fixed common grid built from MAST s_region footprints."""
from __future__ import annotations

import argparse

import numpy as np
import pytest
from astropy.table import Table

from jwst_stack import grid as grid_mod
from jwst_stack.grid import (
    build_grid_wcs,
    load_grid,
    parse_s_region,
    run_grid,
    save_grid,
    summarize_grid,
    union_footprint,
    write_inputs_cache,
)

ROWS = Table(
    {
        "obs_id": ["jw02731-o001_t017_nircam_clear-f200w", "jw02731-o001_t017_nircam_f444w-f470n"],
        "filters": ["F200W", "F444W;F470N"],
        "s_region": [
            "POLYGON 159.13 -58.64 159.13 -58.62 159.15 -58.62 159.15 -58.64",
            "POLYGON 159.10 -58.70 159.10 -58.60 159.16 -58.60 159.16 -58.70",
        ],
    }
)


def test_parse_s_region_polygon():
    v = parse_s_region("POLYGON 159.1 -58.6 159.2 -58.6 159.2 -58.7")
    assert v is not None
    assert v.shape == (3, 2)
    assert v[0, 0] == pytest.approx(159.1)
    assert v[2, 1] == pytest.approx(-58.7)


@pytest.mark.parametrize(
    "value",
    [None, np.ma.masked, "", "CIRCLE 159.1 -58.6 0.1", "NOTAFOOTPRINT"],
)
def test_parse_s_region_rejects(value):
    assert parse_s_region(value) is None


def test_parse_s_region_rejects_odd_tokens():
    assert parse_s_region("POLYGON 159.1 -58.6 159.2") is None


def test_parse_s_region_rejects_bad_numbers():
    assert parse_s_region("POLYGON abc -58.6 159.2 -58.7 159.1 -58.8") is None


def test_union_footprint_bounding_box():
    footprint = union_footprint(ROWS)
    assert footprint.ra_min == pytest.approx(159.10)
    assert footprint.ra_max == pytest.approx(159.16)
    assert footprint.dec_min == pytest.approx(-58.70)
    assert footprint.dec_max == pytest.approx(-58.60)
    assert footprint.rows_used == 2
    assert footprint.rows_total == 2


def test_union_footprint_counts_unusable_rows():
    rows = ROWS.copy()
    rows["s_region"][0] = np.ma.masked
    footprint = union_footprint(rows)
    assert footprint.rows_used == 1
    assert footprint.rows_total == 2


def test_union_footprint_requires_something():
    rows = Table({"obs_id": ["a"], "filters": ["F200W"], "s_region": [np.ma.masked]})
    with pytest.raises(ValueError):
        union_footprint(rows)


def test_build_grid_wcs_north_up_and_covering():
    footprint = union_footprint(ROWS)
    wcs, shape = build_grid_wcs(footprint, pixel_scale_arcsec=0.031, pad_px=10)
    scale_deg = 0.031 / 3600.0
    assert wcs.wcs.cdelt[0] < 0 and wcs.wcs.cdelt[1] > 0
    ny, nx = shape
    assert nx == int(np.ceil((159.16 - 159.10) / scale_deg)) + 20
    assert ny == int(np.ceil((-58.60 - -58.70) / scale_deg)) + 20

    # The grid must strictly contain the footprint.  Because the RA axis is
    # RA---TAN, RA advances cos(dec) * CDELT1 per pixel, so its RA coverage is
    # even broader than the naive bbox/scale estimate.
    corners = np.array(
        [
            [159.10, -58.70],
            [159.16, -58.70],
            [159.16, -58.60],
            [159.10, -58.60],
        ]
    )
    xs, ys = wcs.all_world2pix(corners[:, 0], corners[:, 1], 0)
    assert np.all(xs >= 1) and np.all(xs < nx)
    assert np.all(ys >= 1) and np.all(ys < ny)


def test_save_load_grid_roundtrip(tmp_path):
    footprint = union_footprint(ROWS)
    wcs, shape = build_grid_wcs(footprint)
    path = save_grid(tmp_path / "grid.fits", wcs, shape, 0.031, footprint, ROWS)
    wcs2, shape2 = load_grid(path)
    assert shape2 == shape
    assert str(wcs2.wcs.ctype[0]) == str(wcs.wcs.ctype[0]) == "RA---TAN"
    assert str(wcs2.wcs.ctype[1]) == str(wcs.wcs.ctype[1]) == "DEC--TAN"
    np.testing.assert_allclose(wcs2.wcs.crval, wcs.wcs.crval, rtol=1e-12)
    with open(path, "rb") as handle:
        assert b"GRDNX" in handle.read()


def test_fetch_obs_rows_uses_cache_without_network(monkeypatch, tmp_path):
    cache = tmp_path / "obs.csv"
    write_inputs_cache(ROWS, cache)

    def fail(*a, **k):
        raise AssertionError("network must not be hit")

    monkeypatch.setattr(grid_mod.download, "query_obs_table", fail)
    rows = grid_mod.fetch_obs_rows(2731, cache, refresh=False)
    assert len(rows) == 2
    assert set(rows["filters"]) == {"F200W", "F444W;F470N"}


def test_fetch_obs_rows_queries_and_caches(monkeypatch, tmp_path):
    cache = tmp_path / "obs.csv"

    monkeypatch.setattr(grid_mod.download, "query_obs_table", lambda **kw: ROWS)
    rows = grid_mod.fetch_obs_rows(2731, cache, refresh=False)
    assert len(rows) == 2
    assert cache.exists()


def test_run_grid_builds_then_never_rebuilds(tmp_path, monkeypatch, capsys):
    args = argparse.Namespace(
        outdir=str(tmp_path / "out"),
        scale=0.031,
        pad=10,
        proposal_id=2731,
        refresh_cache=False,
        show=False,
    )

    monkeypatch.setattr(grid_mod.download, "query_obs_table", lambda **kw: ROWS)

    code = run_grid(args)
    assert code == 0
    grid_path = tmp_path / "out" / "grid.fits"
    assert grid_path.exists()
    assert (tmp_path / "out" / "grid_inputs" / "ngc3324_obs.csv").exists()
    before = grid_path.read_bytes()

    def fail(*a, **k):
        raise AssertionError("rebuild would hit the network")

    monkeypatch.setattr(grid_mod.download, "query_obs_table", fail)
    code = run_grid(args)
    assert code == 0
    assert grid_path.read_bytes() == before
    out = capsys.readouterr().out
    assert "common output grid" in out
    assert "pixel scale" in out


def test_summarize_grid_reports_counts(tmp_path):
    footprint = union_footprint(ROWS)
    wcs, shape = build_grid_wcs(footprint)
    summary = summarize_grid(wcs, shape, footprint, ROWS, 0.031)
    assert f"{shape[0]} x {shape[1]} pixels" in summary
    assert "0.0310 arcsec/px" in summary
    assert "F200W" in summary
    assert "F444W;F470N" in summary
    assert "1 obs" in summary