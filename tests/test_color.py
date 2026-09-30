from __future__ import annotations

import json

import numpy as np
import pytest

from jwst_stack.color import (
    CHANNELS,
    COLOR_REFERENCE,
    background_centred_asinh,
    channel_table,
    default_softening,
    feature_channels,
    format_registration_summary_color,
    mosaic_paths,
    rgb_channels,
    shared_limits,
    solve_cross_filter_registration,
    source_free_level,
    stretch_channels,
    subtract_background,
    write_color_registration_json,
)


def _scene(seed: int = 0, ny: int = 128, nx: int = 128) -> np.ndarray:
    """A flat pedestal plus noise plus a few bright sources."""
    rng = np.random.default_rng(seed)
    img = rng.normal(10.0, 1.0, size=(ny, nx))
    img[10:14, 10:14] += 500.0
    img[80:82, 100:104] += 900.0
    return img


# --------------------------------------------------------------------------
# source_free_level
# --------------------------------------------------------------------------


def test_source_free_level_recovers_a_flat_pedestal():
    img = _scene()
    level = source_free_level(img)
    assert 9.5 < level < 10.5


def test_source_free_level_is_not_dragged_up_by_bright_sources():
    plain = source_free_level(_scene())
    spiky = _scene().copy()
    spiky[40:44, 40:44] += 5000.0
    assert source_free_level(spiky) == pytest.approx(plain, abs=0.05)


def test_source_free_level_is_translation_covariant():
    img = _scene()
    assert source_free_level(img + 7.25) == pytest.approx(source_free_level(img) + 7.25)


def test_source_free_level_ignores_nan_holes():
    img = _scene()
    holed = img.copy()
    holed[:64, :64] = np.nan
    assert source_free_level(holed) == pytest.approx(source_free_level(img), abs=0.1)


def test_source_free_level_of_all_nan_is_nan():
    assert np.isnan(source_free_level(np.full((32, 32), np.nan)))


def test_source_free_level_streams_over_a_memmap(tmp_path):
    img = _scene()
    path = tmp_path / "scene.npy"
    np.save(path, img)
    mm = np.load(path, mmap_mode="r")
    assert source_free_level(mm) == pytest.approx(source_free_level(img))


# --------------------------------------------------------------------------
# subtract_background
# --------------------------------------------------------------------------


def test_subtract_background_preserves_nan_and_dtype():
    img = np.array([[1.0, np.nan], [3.0, 5.0]], dtype=np.float32)
    out = subtract_background(img, 1.0)
    assert out.dtype == np.float32
    assert out[0, 0] == 0.0
    assert np.isnan(out[0, 1])


# --------------------------------------------------------------------------
# shared_limits
# --------------------------------------------------------------------------


def test_shared_limits_pool_every_channel():
    a = np.zeros(1000)
    b = np.full(1000, 100.0)
    lo, hi = shared_limits({"a": a, "b": _sprinkle(b)}, lo_pct=0.0, hi_pct=100.0)
    assert lo == pytest.approx(0.0)
    assert hi == pytest.approx(100.0)


def _sprinkle(v):
    out = v.copy()
    out[0] = 0.0
    out[1] = 100.0
    return out


def test_shared_limits_empty_is_nan():
    lo, hi = shared_limits({})
    assert np.isnan(lo) and np.isnan(hi)


# --------------------------------------------------------------------------
# background_centred_asinh
# --------------------------------------------------------------------------


def test_background_centred_asinh_maps_lo_and_hi():
    values = np.array([-5.0, 0.0, 5.0])
    out = background_centred_asinh(values, background=0.0, a=1.0, lo=-5.0, hi=5.0)
    assert out[0] == pytest.approx(0.0)
    assert out[2] == pytest.approx(1.0)
    assert 0.0 < out[1] < 1.0


def test_background_centred_asinh_keeps_nan_and_clips():
    values = np.array([np.nan, -100.0, 100.0])
    out = background_centred_asinh(values, background=0.0, a=1.0, lo=-1.0, hi=1.0)
    assert np.isnan(out[0])
    assert out[1] == 0.0
    assert out[2] == 1.0


def test_background_centred_asinh_rejects_bad_parameters():
    with pytest.raises(ValueError):
        background_centred_asinh(np.zeros(4), background=0.0, a=0.0, lo=-1.0, hi=1.0)
    with pytest.raises(ValueError):
        background_centred_asinh(np.zeros(4), background=0.0, a=1.0, lo=2.0, hi=1.0)


def test_default_softening_scales_with_range():
    assert default_softening(-10.0, 10.0, 0.1) == pytest.approx(2.0)


# --------------------------------------------------------------------------
# stretch_channels policy
# --------------------------------------------------------------------------


def test_stretch_channels_subtracts_each_channel_own_background():
    a = _scene(seed=1)
    b = _scene(seed=2) + 5.0
    result = stretch_channels({"a": a, "b": b}, tile_px=64)
    assert result.backgrounds["b"] - result.backgrounds["a"] == pytest.approx(5.0, abs=0.05)


def test_stretch_channels_does_not_mutate_inputs():
    a = _scene()
    before = a.copy()
    stretch_channels({"a": a}, tile_px=64)
    assert np.array_equal(a, before)


def test_shared_stretch_keeps_a_brighter_channel_brighter():
    base = np.clip(_scene(), 0.0, None)
    result = stretch_channels({"dim": base, "bright": base * 4.0}, tile_px=64)
    assert result.arrays["bright"].mean() > result.arrays["dim"].mean()


def test_additive_pedestal_does_not_change_colour():
    """The F187N +4.8% / F335M +4.5% guarantee.

    A pedestal of any size added to a channel is removed with that channel's own
    measured level, so the displayed colour is unchanged.  The i2d-relative
    percentages are therefore a calibration note, not a colour term.
    """
    green = _scene(seed=3)
    blue = _scene(seed=4)
    red = _scene(seed=5)

    plain = stretch_channels({"g": green, "b": blue, "r": red}, tile_px=64)

    f187n_like = 0.7185
    f335m_like = 0.2517
    pedestalled = stretch_channels(
        {
            "g": green + f187n_like,
            "b": blue + f335m_like,
            "r": red + 0.0431,
        },
        tile_px=64,
    )

    for name in ("g", "b", "r"):
        np.testing.assert_allclose(
            plain.arrays[name], pedestalled.arrays[name], rtol=1e-9, atol=1e-9
        )
    assert pedestalled.backgrounds["g"] == pytest.approx(
        plain.backgrounds["g"] + f187n_like
    )
    assert pedestalled.backgrounds["b"] == pytest.approx(
        plain.backgrounds["b"] + f335m_like
    )
    assert pedestalled.lo == pytest.approx(plain.lo)
    assert pedestalled.hi == pytest.approx(plain.hi)
    assert pedestalled.a == pytest.approx(plain.a)


def test_stretch_result_as_dict_round_trips_parameters():
    result = stretch_channels({"a": _scene()}, tile_px=64)
    payload = result.as_dict()
    assert set(payload["backgrounds"]) == {"a"}
    assert payload["hi"] > payload["lo"]
    assert payload["softening"] > 0.0


# --------------------------------------------------------------------------
# the channel table
# --------------------------------------------------------------------------


def test_there_are_six_channels_three_in_rgb_and_three_features():
    assert len(CHANNELS) == 6
    assert len({c.name for c in CHANNELS}) == 6
    assert {c.role for c in rgb_channels()} == {"R", "G", "B"}
    assert len(feature_channels()) == 3
    assert all(c.is_feature for c in feature_channels())


def test_the_rgb_assignment_is_the_approved_one():
    roles = {c.role: c.bandpass for c in rgb_channels()}
    assert roles == {"R": "F444W", "G": "F200W", "B": "F090W"}


def test_f200w_is_the_registration_reference_and_is_in_the_rgb():
    assert COLOR_REFERENCE == "F200W"
    reference = channel_table(".")[COLOR_REFERENCE]
    assert reference.role == "G"
    assert not reference.is_feature


def test_the_narrowbands_are_features_and_never_in_the_rgb():
    names = {c.name for c in feature_channels()}
    assert names == {"F187N", "F335M", "F444W;F470N"}
    assert "F187N" not in {c.name for c in rgb_channels()}


def test_f470n_is_separable_from_the_clear_f444w_run():
    """Both read FILTER=F444W; the pupil is what keeps them apart."""
    table = channel_table(".")
    clear = table["F444W"]
    narrow = table["F444W;F470N"]
    assert clear.filter_name == narrow.filter_name == "F444W"
    assert clear.pupil == "CLEAR" and narrow.pupil == "F470N"
    assert clear.mosaic != narrow.mosaic
    assert clear.mosaic.startswith("f444w_")
    assert narrow.mosaic.startswith("f444w_f470n_")


def test_bandpass_names_the_pupil_only_when_it_is_not_clear():
    table = channel_table(".")
    assert table["F090W"].bandpass == "F090W"
    assert table["F444W"].bandpass == "F444W"
    assert table["F444W;F470N"].bandpass == "F444W;F470N"


def test_the_colour_build_consumes_the_fine_grid_mosaics_for_the_long_wave_set():
    """F335M/F444W/F470N have two mosaics each; the colour grid is the fine one."""
    table = channel_table(".")
    for name in ("F335M", "F444W", "F444W;F470N"):
        assert table[name].mosaic.endswith("_0031grid.fits"), name


def test_mosaic_paths_resolves_every_channel(tmp_path):
    paths = mosaic_paths(tmp_path)
    assert set(paths) == {c.name for c in CHANNELS}
    for name, path in paths.items():
        assert path.endswith(channel_table(tmp_path)[name].mosaic)


def test_mosaic_paths_rejects_an_unknown_channel(tmp_path):
    with pytest.raises(KeyError, match="F999W"):
        mosaic_paths(tmp_path, ["F090W", "F999W"])


def test_validation_json_prefers_the_fine_grid_variant(tmp_path):
    (tmp_path / "f335m_all_detectors_mosaic_vs_x.json").write_text("{}")
    (tmp_path / "f335m_all_detectors_mosaic_0031grid_vs_x.json").write_text("{}")
    found = channel_table(tmp_path)["F335M"].validation_json(tmp_path)
    assert "0031grid" in found.name


def test_validation_json_of_an_unrun_channel_says_so(tmp_path):
    found = channel_table(tmp_path)["F090W"].validation_json(tmp_path)
    assert "missing" in found.name


# --------------------------------------------------------------------------
# cross-filter registration
# --------------------------------------------------------------------------


def _catalog(n=500, dx=0.0, dy=0.0, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.uniform(10, 240, n)
    y = rng.uniform(10, 240, n)
    return x + dx, y + dy


def test_cross_filter_registration_recovers_a_known_shift():
    ref = _catalog(seed=1)
    moved = _catalog(seed=1, dx=0.31, dy=-0.18)
    solutions = solve_cross_filter_registration(
        {"F200W": ref, "F090W": moved},
        reference="F200W",
        match_radius_px=1.0,
        grid_scale_arcsec=0.031,
    )
    by_name = {s.channel: s for s in solutions}
    assert by_name["F090W"].dx == pytest.approx(0.31, abs=0.01)
    assert by_name["F090W"].dy == pytest.approx(-0.18, abs=0.01)
    assert by_name["F090W"].n_matched == ref[0].size


def test_cross_filter_registration_marks_the_reference_as_zero():
    ref = _catalog(seed=2)
    solutions = solve_cross_filter_registration({"F200W": ref}, reference="F200W")
    assert len(solutions) == 1
    assert solutions[0].is_reference
    assert solutions[0].dx == 0.0 and solutions[0].dy == 0.0
    assert solutions[0].n_detected == ref[0].size


def test_cross_filter_registration_converts_with_the_grid_scale():
    ref = _catalog(seed=3)
    moved = _catalog(seed=3, dx=0.2, dy=0.0)
    solutions = solve_cross_filter_registration(
        {"F200W": ref, "F090W": moved}, reference="F200W",
        match_radius_px=1.0, grid_scale_arcsec=0.031,
    )
    assert solutions[1].median_arcsec == pytest.approx(0.2 * 0.031, rel=1e-6)


def test_a_sub_threshold_shift_is_recorded_but_not_corrected():
    """The Phase 2 finding, as an assertion: 0.04 px is not a correction."""
    ref = _catalog(seed=4)
    moved = _catalog(seed=4, dx=0.0381, dy=0.0017)
    solutions = solve_cross_filter_registration(
        {"F200W": ref, "F090W": moved}, reference="F200W",
        match_radius_px=1.0, grid_scale_arcsec=0.031,
        negligible_px=0.25,
    )
    payload = solutions[1].as_dict()
    assert payload["dx_px"] == pytest.approx(0.0381, abs=0.005)
    assert payload["shift_negligible"] is True
    assert payload["shift_applied"] is False


def test_a_large_shift_is_flagged_as_needing_a_correction():
    ref = _catalog(seed=5)
    moved = _catalog(seed=5, dx=1.4, dy=0.0)
    solutions = solve_cross_filter_registration(
        {"F200W": ref, "F090W": moved}, reference="F200W",
        match_radius_px=3.0, grid_scale_arcsec=0.031, negligible_px=0.25,
    )
    assert solutions[1].negligible is False
    assert solutions[1].as_dict()["shift_applied"] is False, "reported, never applied here"


def test_cross_filter_registration_requires_the_reference():
    with pytest.raises(KeyError, match="F200W"):
        solve_cross_filter_registration({"F090W": _catalog()}, reference="F200W")


def test_cross_filter_registration_reports_the_fraction_within_half_a_pixel():
    ref = _catalog(seed=6)
    moved = _catalog(seed=6, dx=0.1, dy=0.0)
    exact = solve_cross_filter_registration(
        {"F200W": ref, "F090W": moved}, reference="F200W", match_radius_px=1.0
    )
    assert exact[1].within_half_px == pytest.approx(1.0, abs=0.05)

    jittered = _catalog(seed=7)
    noisy = solve_cross_filter_registration(
        {"F200W": ref, "F090W": jittered}, reference="F200W", match_radius_px=1.0
    )
    assert noisy[1].within_half_px < exact[1].within_half_px


def test_registration_json_records_the_basis_and_why_nothing_shifted(tmp_path):
    ref = _catalog(seed=8)
    moved = _catalog(seed=8, dx=0.04, dy=0.0)
    solutions = solve_cross_filter_registration(
        {"F200W": ref, "F090W": moved}, reference="F200W",
        match_radius_px=1.0, grid_scale_arcsec=0.031,
    )
    out = write_color_registration_json(
        tmp_path / "color_registration_f200w.json",
        solutions,
        grid_scale_arcsec=0.031,
        grid_path="out/grid.fits",
    )
    payload = json.loads(out.read_text(encoding="utf-8"))["cross_filter_registration"]
    assert payload["reference"] == "F200W"
    assert payload["pixel_basis"] == "grid_px"
    assert payload["correction_applied"] is False
    assert "below the threshold" in payload["correction_reason"]
    assert set(payload["channels"]) == {"F200W", "F090W"}


def test_registration_json_reports_the_reference_scale_not_its_zero_shift(tmp_path):
    """The reference's shift is zero by definition; its scale is the grid's.

    Reading the field off the reference solution would write 0.0 into a field
    that claims to be an angular pixel scale - and 0.0 is not a plausible
    scale, so a consumer would have no way to notice.
    """
    ref = _catalog(seed=8)
    solutions = solve_cross_filter_registration(
        {"F200W": ref, "F090W": _catalog(seed=8, dx=0.04, dy=0.0)},
        reference="F200W", match_radius_px=1.0, grid_scale_arcsec=0.031,
    )
    assert solutions[0].channel == "F200W"
    assert solutions[0].median_arcsec == 0.0
    out = write_color_registration_json(
        tmp_path / "reg.json", solutions, grid_scale_arcsec=0.031
    )
    payload = json.loads(out.read_text(encoding="utf-8"))["cross_filter_registration"]
    assert payload["reference_scale_arcsec_per_px"] == pytest.approx(0.031)


def test_registration_json_warns_when_a_channel_needs_shifting(tmp_path):
    ref = _catalog(seed=9)
    moved = _catalog(seed=9, dx=0.9, dy=0.0)
    solutions = solve_cross_filter_registration(
        {"F200W": ref, "F090W": moved}, reference="F200W",
        match_radius_px=2.0, grid_scale_arcsec=0.031,
    )
    out = write_color_registration_json(tmp_path / "reg.json", solutions)
    payload = json.loads(out.read_text(encoding="utf-8"))["cross_filter_registration"]
    assert payload["correction_applied"] is True
    assert "need shifting" in payload["correction_reason"]


def test_registration_summary_labels_the_basis():
    ref = _catalog(seed=10)
    moved = _catalog(seed=10, dx=0.05, dy=0.0)
    solutions = solve_cross_filter_registration(
        {"F200W": ref, "F090W": moved}, reference="F200W", match_radius_px=1.0
    )
    text = format_registration_summary_color(solutions, 0.031)
    assert "grid px" in text
    assert "reference" in text
    assert "F090W" in text


# --------------------------------------------------------------------------
# the CLI surface
# --------------------------------------------------------------------------


def test_color_subcommands_parse():
    from jwst_stack.cli import build_parser

    parser = build_parser()
    psf = parser.parse_args(["color-psf", "--channels", "F090W", "--tiles", "0"])
    assert psf.channels == ["F090W"]
    assert psf.tiles == 0, "0 means every covered tile"
    assert psf.out.endswith("color_channels.json")

    reg = parser.parse_args(
        ["color-register", "--reference", "F200W", "--match-radius-px", "4"]
    )
    assert reg.reference == "F200W"
    assert reg.match_radius_px == 4.0
    assert reg.out.endswith("color_registration_f200w.json")


def test_color_psf_defaults_to_refining():
    from jwst_stack.cli import build_parser

    assert build_parser().parse_args(["color-psf"]).no_refine is False
    assert build_parser().parse_args(["color-psf", "--no-refine"]).no_refine is True
