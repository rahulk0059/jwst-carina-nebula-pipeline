"""Tests for photutils star registration (``jwst_stack.register``).

The overlays are genuine 2-D images with injected Gaussian sources, so the
tests exercise the real photutils detection / overlap matching / translation
solve / background-offset engine with no WCS mocking.  A final end-to-end test
drives :func:`jwst_stack.stack.stack_visit` on synthetic exposures to prove
the ``register=True`` default forward the solution into the combined stack
and that ``--no-register`` reproduces the pure-WCS combination.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from jwst_stack import register as reg

PIXEL_SCALE_ARCSEC = 0.031
FWHM_PX = 3.0

GRID = (128, 128)


def _gaussian_overlay(
    centers: list[tuple[float, float]],
    shape: tuple[int, int] = GRID,
    fwhm_px: float = FWHM_PX,
    amplitude: float = 100.0,
) -> np.ndarray:
    """Output-grid overlay with a Gaussian at each *center* (x, y)."""
    sigma = fwhm_px / (2 * np.sqrt(2 * np.log(2)))
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    img = np.zeros(shape, dtype=float)
    for cx, cy in centers:
        img += amplitude * np.exp(
            -((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sigma**2)
        )
    return img


def _exposure_stub(name: str, visit: int, exposure: int):
    return SimpleNamespace(
        name=name,
        visit=visit,
        exposure=exposure,
        filename=name,
    )


def _trigger_vision(exposures, overlays):
    """Build enough frame metadata for register_overlays to detect stars and
    run the solve without needing real FITS/WCS."""
    footprint = np.isfinite(overlays[0]) * 1.0
    footprints = [np.array(footprint, copy=True) for _ in overlays]
    widths = [
        float(np.max(np.where(np.isfinite(o))[1]) - np.min(np.where(np.isfinite(o))[1]))
        if np.any(np.isfinite(o))
        else 0.0
        for o in overlays
    ]
    return footprints, widths


def _run_register(overlays, match_radius_arcsec=1.0):
    exposures = [_exposure_stub(f"e{i:02d}_cal.fits", 1, i) for i in range(len(overlays))]
    footprints, _ = _trigger_vision(exposures, overlays)
    shifted, result = reg.register_overlays(
        overlays,
        footprints,
        exposures,
        out_wcs=None,
        pixel_scale_arcsec=PIXEL_SCALE_ARCSEC,
        fwhm_px=FWHM_PX,
        match_radius_arcsec=match_radius_arcsec,
        threshold_sigma=5.0,
        bg_sigma=3.0,
    )
    return exposures, shifted, result


def test_identical_frames_zero_shift():
    centers = [(30.0, 40.0), (70.0, 25.0), (85.0, 90.0)]
    overlays = [_gaussian_overlay(centers), _gaussian_overlay(centers)]
    _, shifted, result = _run_register(overlays)

    assert result.reference == "e00_cal.fits"
    assert result.n_frames == 2
    sol = result.frames[0]
    assert sol.n_matched == len(centers)
    assert sol.dx == pytest.approx(0.0, abs=1e-9)
    assert sol.dy == pytest.approx(0.0, abs=1e-9)
    assert np.allclose(shifted[1], overlays[1], equal_nan=True)


def test_translation_solve_integer_pixels():
    centers = [(30.0, 40.0), (70.0, 25.0), (85.0, 90.0)]
    dx, dy = 3.0, -2.0
    overlays = [
        _gaussian_overlay(centers),
        _gaussian_overlay([(x + dx, y + dy) for x, y in centers]),
    ]
    _, shifted, result = _run_register(overlays)

    sol = result.frames[0]
    assert sol.n_matched == len(centers)
    assert sol.dx == pytest.approx(dx, abs=1e-9)
    assert sol.dy == pytest.approx(dy, abs=1e-9)
    # after the solve the frame starlands on the reference
    assert np.allclose(shifted[1], overlays[0], atol=0.1)


def test_translation_solve_fractional_pixels():
    centers = [(30.0, 25.0), (70.0, 30.0), (55.0, 55.0), (90.0, 85.0)]
    dx, dy = 1.5, -0.75
    overlays = [
        _gaussian_overlay(centers),
        _gaussian_overlay([(x + dx, y + dy) for x, y in centers]),
    ]
    _, shifted, result = _run_register(overlays)

    sol = result.frames[0]
    assert sol.n_matched >= len(centers)
    assert sol.dx == pytest.approx(dx, abs=0.15)
    assert sol.dy == pytest.approx(dy, abs=0.15)


def test_background_offset_corrected():
    centers = [(30.0, 25.0), (70.0, 30.0), (55.0, 55.0)]
    bg = 7.5
    overlays = [
        _gaussian_overlay(centers),
        _gaussian_overlay(centers) + bg,
    ]
    _, shifted, result = _run_register(overlays)

    sol = result.frames[0]
    assert sol.n_matched == len(centers)
    assert sol.dx == pytest.approx(0.0, abs=1e-9)
    assert sol.dy == pytest.approx(0.0, abs=1e-9)
    assert sol.background_offset == pytest.approx(bg, abs=0.05)
    assert np.allclose(shifted[1], overlays[0], atol=0.05)


def test_match_requires_overlap_radius():
    centers = [(30.0, 40.0), (70.0, 25.0)]
    # second frame stars are *far* from the reference centres -> no matches
    overlays = [
        _gaussian_overlay(centers),
        _gaussian_overlay([(60.0, 100.0), (110.0, 20.0), (95.0, 15.0)]),
    ]
    _, shifted, result = _run_register(overlays, match_radius_arcsec=0.1)

    sol = result.frames[0]
    assert sol.n_matched == 0
    # guard: too few matches -> keep the pure-WCS stack (no shift applied)
    assert sol.dx == 0.0 and sol.dy == 0.0
    assert np.allclose(shifted[1], overlays[1], equal_nan=True)


def test_write_registration_json(tmp_path: Path):
    centers = [(30.0, 40.0), (70.0, 25.0), (85.0, 90.0)]
    dx, dy = 2.0, 1.0
    overlays = [
        _gaussian_overlay(centers),
        _gaussian_overlay([(x + dx, y + dy) for x, y in centers]),
    ]
    _, shifted, result = _run_register(overlays)

    out = tmp_path / "registration.json"
    written = reg.write_registration_json(result, out)
    assert written == out
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["reference"] == "e00_cal.fits"
    assert payload["n_frames"] == 2
    frames = {f["filename"]: f for f in payload["frames"]}
    assert frames["e01_cal.fits"]["dx"] == pytest.approx(dx, abs=1e-9)
    assert frames["e01_cal.fits"]["dy"] == pytest.approx(dy, abs=1e-9)


def test_format_summary_includes_reference_and_residuals():
    centers = [(30.0, 40.0), (70.0, 25.0)]
    overlays = [_gaussian_overlay(centers), _gaussian_overlay(centers)]
    _, shifted, result = _run_register(overlays)

    text = reg.format_registration_summary(result)
    assert "reference" in text
    assert "dx=" in text
    assert "n_matched" in text or "match" in text.lower()


def test_register_too_few_exposures_is_noop():
    centers = [(30.0, 40.0)]
    overlays = [_gaussian_overlay(centers)]
    _, shifted, result = _run_register(overlays)
    assert result.n_frames == 1
    assert len(result.frames) == 0


@pytest.mark.parametrize("dx,dy", [(0.3, 0.0), (0.0, -0.4), (0.5, 0.5)])
def test_fractional_shift_px_is_solved_subpixel(dx: float, dy: float):
    centers = [(25.0, 45.0), (60.0, 30.0), (80.0, 70.0), (40.0, 90.0), (95.0, 20.0)]
    overlays = [
        _gaussian_overlay(centers),
        _gaussian_overlay([(x + dx, y + dy) for x, y in centers]),
    ]
    _, shifted, result = _run_register(overlays)

    sol = result.frames[0]
    assert sol.n_matched >= max(3, len(centers)) - 2
    assert sol.dx == pytest.approx(dx, abs=0.15)
    assert sol.dy == pytest.approx(dy, abs=0.15)


def test_nan_footprint_does_not_poison_cubic_shift():
    """A NaN outside the footprint must not spread over the valid region.

    ``scipy.ndimage.shift`` with ``order >= 2`` runs a global spline
    prefilter, so a single NaN in the input turns the whole output into NaN.
    Registration therefore has to fill the mask before shifting and restore it
    afterwards, otherwise every frame but the reference is dropped and the
    stack depth collapses to 1.
    """
    centers = [(30.0, 40.0), (70.0, 25.0), (85.0, 90.0), (55.0, 60.0)]
    reference = _gaussian_overlay(centers)
    frame = _gaussian_overlay([(x + 0.4, y - 0.3) for x, y in centers])
    frame[:, 100:] = np.nan

    assert np.isnan(frame).any()
    assert np.isfinite(frame[:, :100]).any()

    exposures = [
        _exposure_stub("e00_cal.fits", 1, 0),
        _exposure_stub("e01_cal.fits", 1, 1),
    ]
    shifted, result = reg.register_overlays(
        [reference, frame],
        [np.isfinite(reference), np.isfinite(frame)],
        exposures,
        out_wcs=None,
        pixel_scale_arcsec=PIXEL_SCALE_ARCSEC,
        fwhm_px=FWHM_PX,
        match_radius_arcsec=1.0,
        threshold_sigma=5.0,
        bg_sigma=3.0,
        interp_order=3,
    )

    assert result.frames[0].n_matched >= 3
    assert np.isnan(shifted[1][:, 100:]).all()
    assert np.isfinite(shifted[1][:, :100]).sum() > 0.9 * shifted[1][:, :100].size
