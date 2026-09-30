from __future__ import annotations

import json
import math

import numpy as np
import pytest

from jwst_stack.io import read_json_section, update_json_section
from jwst_stack.psf import (
    FWHM_PER_SIGMA,
    fit_star_patch,
    gaussian_2d,
    measure_fwhm,
    measure_channels_fwhm,
    patch_half_width,
    write_fwhm_json,
)


def _star_patch(fwhm_px: float, half: int, amp: float = 500.0, pedestal: float = 7.0, theta=0.3):
    """A star centred in a patch, plus its centre in patch coordinates."""
    size = 2 * half + 1
    yy, xx = np.mgrid[0:size, 0:size]
    xy = np.stack([xx.ravel(), yy.ravel()])
    model = gaussian_2d(
        xy, pedestal, amp, half, half, fwhm_px / FWHM_PER_SIGMA,
        fwhm_px / FWHM_PER_SIGMA, theta,
    )
    return model.reshape(size, size), half, half


# --------------------------------------------------------------------------
# the model and the patch size
# --------------------------------------------------------------------------


def test_gaussian_2d_peaks_at_its_centre():
    half = 6
    patch, cx, cy = _star_patch(4.0, half)
    assert patch[cy, cx] == pytest.approx(patch.max())


def test_gaussian_2d_halves_at_the_fwhm_radius():
    """A width definition that is off by a constant would be invisible
    everywhere else in the pipeline, so it is pinned here directly."""
    sigma = 1.7
    half = 20
    size = 2 * half + 1
    yy, xx = np.mgrid[0:size, 0:size]
    xy = np.stack([xx.ravel(), yy.ravel()])
    model = gaussian_2d(xy, 0.0, 1.0, half, half, sigma, sigma, 0.0).reshape(size, size)
    peak = model[half, half]
    half_max_radius = round(sigma * FWHM_PER_SIGMA / 2)
    assert model[half, half + half_max_radius] == pytest.approx(peak / 2, rel=0.02)


def test_patch_half_width_grows_with_the_psf():
    assert patch_half_width(2.0) == 6
    assert patch_half_width(3.0) == 6, "narrow PSFs hit the floor"
    assert patch_half_width(5.5) == 10
    # 4 sigma at 5.5 px FWHM, so the patch contains the wings.
    assert patch_half_width(5.5) >= 4 * 5.5 / FWHM_PER_SIGMA


# --------------------------------------------------------------------------
# the per-star fit
# --------------------------------------------------------------------------


@pytest.mark.parametrize("fwhm", [2.0, 3.5, 5.5])
def test_fit_star_patch_recovers_a_known_fwhm(fwhm):
    patch, cx, cy = _star_patch(fwhm, patch_half_width(fwhm))
    fit = fit_star_patch(patch, cx, cy, fwhm_hint_px=3.0)
    assert fit is not None
    assert fit.fwhm_px == pytest.approx(fwhm, rel=0.01)


def test_fit_star_patch_recovers_a_pedestal():
    """The pedestal is fitted, never subtracted from an external level."""
    patch, cx, cy = _star_patch(3.0, patch_half_width(3.0), pedestal=12.5)
    fit = fit_star_patch(patch, cx, cy)
    assert fit.pedestal == pytest.approx(12.5, abs=0.1)


def test_fit_star_patch_is_invariant_to_an_additive_offset():
    """A pedestal difference between channels cannot move the fitted width.

    This is the same invariance the colour policy rests on, tested at the level
    where it is actually implemented: a constant added to the pixel data changes
    only the fitted pedestal.
    """
    patch, cx, cy = _star_patch(3.0, patch_half_width(3.0))
    a = fit_star_patch(patch, cx, cy)
    b = fit_star_patch(patch + 11.0, cx, cy)
    assert a.fwhm_px == pytest.approx(b.fwhm_px, rel=0.01)
    assert b.pedestal == pytest.approx(a.pedestal + 11.0, abs=0.1)


def test_fit_star_patch_reports_both_axes_of_an_elliptical_psf():
    size = 2 * patch_half_width(5.5) + 1
    yy, xx = np.mgrid[0:size, 0:size]
    xy = np.stack([xx.ravel(), yy.ravel()])
    sig = 5.5 / FWHM_PER_SIGMA
    patch = gaussian_2d(xy, 0.0, 500.0, size / 2, size / 2, sig, sig * 0.5, 0.0)
    fit = fit_star_patch(patch.reshape(size, size), size / 2, size / 2)
    assert fit.fwhm_major_px == pytest.approx(5.5, rel=0.02)
    assert fit.fwhm_minor_px == pytest.approx(2.75, rel=0.03)
    assert fit.fwhm_px == pytest.approx(math.sqrt(5.5 * 2.75), rel=0.03)
    assert fit.axis_ratio == pytest.approx(2.0, rel=0.05)


def test_fit_star_patch_returns_none_on_a_flat_patch():
    flat = np.full((13, 13), 5.0)
    assert fit_star_patch(flat, 6, 6) is None


def test_fit_star_patch_returns_none_on_an_unphysical_width():
    """A fit that lands outside the plausible range is dropped, not averaged in."""
    half = 20
    patch, cx, cy = _star_patch(4.0, half)
    patch[cy - 1 : cy + 2, cx - 1 : cx + 2] = 1e6  # a saturated 3x3 core
    fit = fit_star_patch(patch, cx, cy)
    assert fit is None or fit.fwhm_px < 30.0


def test_fit_star_patch_rejects_a_patch_with_nan():
    patch, cx, cy = _star_patch(3.0, 8)
    patch[0, 0] = np.nan
    assert fit_star_patch(patch, cx, cy) is None


# --------------------------------------------------------------------------
# the channel aggregate
# --------------------------------------------------------------------------


def _star_grid(shape=(256, 256), fwhm=3.0, amp=900.0, pedestal=10.0, seed=0, step=32):
    rng = np.random.default_rng(seed)
    img = np.full(shape, pedestal) + rng.normal(0.0, 0.5, shape)
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    sig = fwhm / FWHM_PER_SIGMA
    for cy in range(step, shape[0] - step, step):
        for cx in range(step, shape[1] - step, step):
            img += amp * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sig**2))
    return img


def test_measure_fwhm_recovers_a_planted_width():
    img = _star_grid(fwhm=3.0)
    result = measure_fwhm(
        img, None, img.shape, channel="T", grid_scale_arcsec=0.031,
        n_tiles=4, tile_px=128, target_stars=40, separation_px=12.0,
        fwhm_hint_px=3.0,
    )
    assert result.measured
    assert result.fwhm_px == pytest.approx(3.0, rel=0.03)
    assert result.fwhm_arcsec == pytest.approx(0.031 * result.fwhm_px, rel=1e-9)


def test_measure_fwhm_is_unmoved_by_an_additive_pedestal():
    base = _star_grid(fwhm=3.0)
    kwargs = dict(
        n_tiles=4, tile_px=128, target_stars=40, fwhm_hint_px=3.0,
        grid_scale_arcsec=0.031,
    )
    a = measure_fwhm(base, None, base.shape, **kwargs)
    b = measure_fwhm(base + 20.0, None, base.shape, **kwargs)
    assert a.fwhm_px == pytest.approx(b.fwhm_px, rel=0.01)


def test_measure_fwhm_of_a_blank_field_is_not_measured():
    blank = np.full((128, 128), 5.0)
    result = measure_fwhm(blank, None, blank.shape, n_tiles=2, tile_px=64)
    assert not result.measured
    assert np.isnan(result.fwhm_px)
    assert result.n_used == 0


def test_measure_fwhm_of_a_sparse_field_reports_not_measured():
    """Below the star floor the channel is reported as unmeasured, not as zero."""
    img = np.full((256, 256), 5.0)
    img[100:110, 100:110] += 500.0
    result = measure_fwhm(img, None, img.shape, n_tiles=4, tile_px=128, target_stars=40)
    assert not result.measured


def test_measure_fwhm_counts_the_tiles_it_used():
    img = _star_grid(shape=(256, 256))
    result = measure_fwhm(img, None, img.shape, n_tiles=3, tile_px=128, fwhm_hint_px=3.0)
    assert result.n_tiles == 3


# --------------------------------------------------------------------------
# the two passes, and the record they write
# --------------------------------------------------------------------------


def test_measure_channels_fwhm_returns_both_passes(tmp_path):
    from astropy.io import fits

    from .conftest import make_wcs
    from tests.test_starcat import _write_mosaic

    wcs = make_wcs(crpix=(64.5, 64.5))
    paths = {}
    for name, fwhm in (("A", 2.5), ("B", 4.0)):
        p = tmp_path / f"{name}.fits"
        _write_mosaic(p, _star_grid(fwhm=fwhm, step=24), wcs=wcs)
        paths[name] = str(p)
    refined, first = measure_channels_fwhm(
        paths, grid_scale_arcsec=0.031, n_tiles=2, tile_px=128,
        target_stars=30, fwhm_hint_px=3.0,
    )
    assert set(refined) == set(first) == {"A", "B"}
    assert refined["A"].fwhm_px == pytest.approx(2.5, rel=0.05)
    assert refined["B"].fwhm_px == pytest.approx(4.0, rel=0.05)
    assert first["A"].detect_fwhm_px == 3.0
    assert refined["A"].detect_fwhm_px == pytest.approx(first["A"].fwhm_px)


def test_write_fwhm_json_preserves_other_sections(tmp_path):
    """The backgrounds and stretch land in the same file, so the writers must
    not truncate each other."""
    record = tmp_path / "color_channels.json"
    update_json_section(record, "backgrounds", {"F090W": 1.2134})
    img = _star_grid(fwhm=3.0)
    result = measure_fwhm(
        img, None, img.shape, channel="F090W", grid_scale_arcsec=0.031,
        n_tiles=4, tile_px=128, target_stars=40, fwhm_hint_px=3.0,
    )
    write_fwhm_json(record, {"F090W": result}, grid_path="out/grid.fits",
                    grid_shape=(256, 256))
    payload = json.loads(record.read_text(encoding="utf-8"))
    assert payload["backgrounds"] == {"F090W": 1.2134}
    assert payload["fwhm"]["channels"]["F090W"]["fwhm_px"] == pytest.approx(3.0, rel=0.05)
    assert payload["fwhm"]["channels"]["F090W"]["pixel_basis"] == "grid_px"
    assert payload["fwhm"]["grid"]["shape"] == [256, 256]


def test_write_fwhm_json_records_the_first_pass_when_given(tmp_path):
    img = _star_grid(fwhm=3.0)
    kwargs = dict(n_tiles=4, tile_px=128, target_stars=40, grid_scale_arcsec=0.031)
    first = measure_fwhm(img, None, img.shape, channel="X", fwhm_hint_px=3.0, **kwargs)
    second = measure_fwhm(img, None, img.shape, channel="X", fwhm_hint_px=3.1, **kwargs)
    record = tmp_path / "c.json"
    write_fwhm_json(record, {"X": second}, first_pass={"X": first})
    payload = json.loads(record.read_text(encoding="utf-8"))
    assert payload["fwhm"]["pass1_detect_hint_3px"]["X"]["fwhm_px"] == pytest.approx(3.0, rel=0.05)
    assert payload["fwhm"]["channels"]["X"]["detect_fwhm_px"] == pytest.approx(3.1)

    # Without a first pass the key is absent rather than present-and-empty, so
    # a reader cannot mistake a single-pass run for a refined one.
    single = tmp_path / "single.json"
    write_fwhm_json(single, {"X": second})
    assert "pass1_detect_hint_3px" not in read_json_section(single, "fwhm")


def test_update_json_section_replaces_a_section_wholesale(tmp_path):
    record = tmp_path / "x.json"
    update_json_section(record, "s", {"a": 1, "b": 2})
    update_json_section(record, "s", {"a": 9})
    assert json.loads(record.read_text(encoding="utf-8"))["s"] == {"a": 9}


def test_read_json_section_of_a_missing_file_is_none(tmp_path):
    assert read_json_section(tmp_path / "nope.json", "s") is None


def test_read_json_section_of_corrupt_json_is_none(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", encoding="utf-8")
    assert read_json_section(bad, "s") is None
