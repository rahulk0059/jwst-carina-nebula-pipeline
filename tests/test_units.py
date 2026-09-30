"""Tests for :mod:`jwst_stack.units`.

Two bugs motivate this module, and both are pinned here as regressions:

* degrees were multiplied by 1e6 instead of 3600, turning a 0.24 arcsec
  disagreement into a bogus 67.7 arcsec;
* an i2d-pixel offset was converted with the 0.031 *grid* scale instead of the
  filter's own *i2d* scale, understating every long-wave cross-check by 50.7%
  and manufacturing a false "the two grids agree" result.

Neither is caught by a round-trip test, so the tests below assert the absolute
values against the numbers recorded in AGENTS.md, and assert that
:func:`px_to_arcsec` cannot be called without a basis.
"""

from __future__ import annotations

import numpy as np
import pytest

from jwst_stack.units import (
    ARCSEC_PER_DEG,
    GRID_PX,
    I2D_PX,
    PIXEL_BASES,
    arcsec_to_deg,
    deg_to_arcsec,
    format_pair_px,
    format_px,
    px_label,
    px_to_arcsec,
)

#: Measured i2d native scales.  Deliberately the real, unrounded values: the
#: whole failure mode is that these are not 0.031 and not 0.0629.
I2D_SCALE_SW = 0.031227
I2D_SCALE_F335M = 0.062904
I2D_SCALE_F444W = 0.062908
I2D_SCALE_F470N = 0.062936

#: The two output grids.
GRID_SCALE_FINE = 0.031
GRID_SCALE_LONGWAVE = 0.0629


def test_degrees_to_arcsec_is_times_3600_not_1e6():
    assert deg_to_arcsec(1.0) == pytest.approx(3600.0)
    assert deg_to_arcsec(1.0) == pytest.approx(1.0 * ARCSEC_PER_DEG)
    # The recorded 0.24 arcsec s_region disagreement between F335M and F470N.
    assert deg_to_arcsec(1e-6) == pytest.approx(0.0036)
    assert deg_to_arcsec(6.666666666666667e-8) == pytest.approx(0.00024)
    # And it is emphatically not microarcsec.
    assert deg_to_arcsec(1e-6) != pytest.approx(1.0)
    assert deg_to_arcsec(1e-6) < 1.0


def test_deg_to_arcsec_is_not_a_microarcsec_conversion():
    """The literal bug: x1e6 instead of x3600, a factor of 277.8."""
    assert deg_to_arcsec(1e-6) == pytest.approx(0.0036)
    assert deg_to_arcsec(1e-6) * 1e6 != pytest.approx(0.0036)


def test_arcsec_deg_round_trip():
    for v in (0.0036, 0.00024, 0.0629, 1.0, 3600.0):
        assert arcsec_to_deg(deg_to_arcsec(v)) == pytest.approx(v)
    assert arcsec_to_deg(1.0) == pytest.approx(1.0 / 3600.0)


def test_px_to_arcsec_requires_a_basis():
    """``basis`` has no default on purpose, so this cannot regress silently."""
    with pytest.raises(TypeError):
        px_to_arcsec(0.0818, 0.062936)  # type: ignore[call-arg]


def test_px_to_arcsec_rejects_an_unknown_basis():
    with pytest.raises(ValueError, match="unknown pixel basis"):
        px_to_arcsec(1.0, 0.031, "sky_px")
    with pytest.raises(ValueError, match="unknown pixel basis"):
        px_to_arcsec(1.0, 0.031, "I2D")  # case matters
    with pytest.raises(ValueError, match="unknown pixel basis"):
        px_to_arcsec(1.0, 0.031, "")


def test_px_to_arcsec_rejects_a_nonsense_scale():
    for bad in (0.0, -0.031, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="pixel scale"):
            px_to_arcsec(1.0, bad, GRID_PX)


def test_px_to_arcsec_known_bases():
    assert px_to_arcsec(0.0818, I2D_SCALE_F470N, I2D_PX) == pytest.approx(0.00515, abs=5e-6)
    assert px_to_arcsec(0.1801, I2D_SCALE_F335M, I2D_PX) == pytest.approx(0.01133, abs=5e-6)
    assert px_to_arcsec(0.1923, I2D_SCALE_F444W, I2D_PX) == pytest.approx(0.01210, abs=5e-6)


def test_px_to_arcsec_vectorises():
    got = px_to_arcsec(np.array([0.0818, 0.1801, 0.1923]), I2D_SCALE_F335M, I2D_PX)
    assert got.shape == (3,)
    assert got[0] == pytest.approx(0.0818 * I2D_SCALE_F335M)


def test_wrong_basis_is_a_50_percent_error_on_the_long_wave_cross_checks():
    """The measured magnitude of the documented bug.

    An i2d-pixel offset converted with the 0.031 grid scale instead of the
    ~0.0629 i2d scale comes out 50.7% too small for the long-wave filters.
    """
    offset_px = 0.1492  # F470N 0.031-grid cross-check, i2d pixels
    correct = px_to_arcsec(offset_px, I2D_SCALE_F470N, I2D_PX)
    wrong = px_to_arcsec(offset_px, GRID_SCALE_FINE, GRID_PX)
    assert correct == pytest.approx(0.00939, abs=5e-6)
    assert wrong / correct == pytest.approx(GRID_SCALE_FINE / I2D_SCALE_F470N)
    assert wrong == pytest.approx(correct * 0.4926, rel=1e-3)
    # The error is the 2.029x grid ratio, i.e. ~50.7% understated.
    assert I2D_SCALE_F470N / GRID_SCALE_FINE == pytest.approx(2.0290, rel=1e-3)


def test_the_four_i2d_scales_are_all_distinct():
    scales = [I2D_SCALE_SW, I2D_SCALE_F335M, I2D_SCALE_F444W, I2D_SCALE_F470N]
    assert len(set(scales)) == 4
    # Short-wave and long-wave i2ds differ by ~2.03x, which is why the bases
    # must be named rather than assumed.
    assert I2D_SCALE_F470N / I2D_SCALE_SW == pytest.approx(2.0146, rel=1e-3)


def test_grid_basis_uses_the_grid_scale():
    assert px_to_arcsec(0.5, GRID_SCALE_FINE, GRID_PX) == pytest.approx(0.0155)
    assert px_to_arcsec(0.5, GRID_SCALE_LONGWAVE, GRID_PX) == pytest.approx(0.03145)


def test_px_label():
    assert px_label(I2D_PX) == "i2d px"
    assert px_label(GRID_PX) == "grid px"
    with pytest.raises(ValueError, match="unknown pixel basis"):
        px_label("nope")


def test_pixel_bases_is_the_closed_registry():
    assert PIXEL_BASES == (I2D_PX, GRID_PX)


def test_format_px_names_basis_and_scale():
    got = format_px(0.0818, I2D_SCALE_F470N, I2D_PX)
    assert "i2d_px" in got
    assert "0.062936" in got
    assert "0.00515" in got
    assert "arcsec" in got


def test_format_px_handles_nan():
    assert format_px(float("nan"), I2D_SCALE_F470N, I2D_PX) == "n/a"
    assert format_px(float("inf"), I2D_SCALE_F470N, I2D_PX) == "n/a"


def test_format_px_rejects_an_unknown_basis():
    with pytest.raises(ValueError, match="unknown pixel basis"):
        format_px(0.1, 0.031, "mosaic_px")


def test_format_pair_px():
    got = format_pair_px(0.1492, -0.0012, I2D_SCALE_F470N, I2D_PX)
    assert "dx +0.1492" in got and "dy -0.0012" in got
    assert "grid_px" not in got
    assert "i2d_px" in got
    assert format_pair_px(float("nan"), 0.0, 0.031, GRID_PX) == "n/a"


def test_format_pair_px_signs_survive():
    got = format_pair_px(-0.0333, 0.0059, GRID_SCALE_FINE, GRID_PX)
    assert "dx -0.0333" in got
    assert "dy +0.0059" in got
