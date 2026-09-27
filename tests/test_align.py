from __future__ import annotations

import numpy as np
import pytest

from jwst_stack.align import build_output_wcs, reproject_to_grid
from tests.conftest import (
    centroid_above,
    gaussian_frame,
    make_exposures,
    make_wcs,
    shift_crval_to_match,
)


def test_reproject_identity_single_frame():
    wcs = make_wcs()
    data = gaussian_frame(wcs, peak_px=(30.0, 30.0))
    out_wcs, shape = build_output_wcs(make_exposures(wcs))
    assert shape == out_wcs.array_shape
    reproj, footprint = reproject_to_grid(data, wcs, out_wcs, method="interp")
    centroid = centroid_above(reproj, threshold=0.5)
    assert centroid is not None
    # Single-frame grid keeps the 32 px margin, so the peak lands at 30 + 32.
    assert centroid[0] == pytest.approx(62.0, abs=0.5)
    assert centroid[1] == pytest.approx(62.0, abs=0.5)


def test_two_frames_align_on_common_grid():
    wcs_a = make_wcs()
    # Frame B moves the same source from (30, 30) to (33, 35) pixels.
    wcs_b = shift_crval_to_match(wcs_a, (30.0, 30.0), (33.0, 35.0))

    data_a = gaussian_frame(wcs_a, peak_px=(30.0, 30.0))
    data_b = gaussian_frame(wcs_b, peak_px=(33.0, 35.0))

    out_wcs, _ = build_output_wcs(make_exposures(wcs_a, wcs_b))
    ra, _ = reproject_to_grid(data_a, wcs_a, out_wcs, method="interp")
    rb, _ = reproject_to_grid(data_b, wcs_b, out_wcs, method="interp")

    ca, cb = centroid_above(ra), centroid_above(rb)
    assert ca is not None and cb is not None
    # Same sky point must land on the same output pixel.
    assert abs(ca[0] - cb[0]) < 0.7
    assert abs(ca[1] - cb[1]) < 0.7


def test_exact_and_interp_agree():
    wcs_a = make_wcs()
    wcs_b = shift_crval_to_match(wcs_a, (30.0, 30.0), (33.0, 35.0))
    data_b = gaussian_frame(wcs_b, peak_px=(33.0, 35.0))
    out_wcs, _ = build_output_wcs(make_exposures(wcs_a, wcs_b))

    r_interp, _ = reproject_to_grid(data_b, wcs_b, out_wcs, method="interp")
    r_exact, _ = reproject_to_grid(data_b, wcs_b, out_wcs, method="exact")

    assert centroid_above(r_interp) is not None
    assert centroid_above(r_exact) is not None
    assert centroid_above(r_interp)[0] == pytest.approx(
        centroid_above(r_exact)[0], abs=0.5
    )
    assert centroid_above(r_interp)[1] == pytest.approx(
        centroid_above(r_exact)[1], abs=0.5
    )