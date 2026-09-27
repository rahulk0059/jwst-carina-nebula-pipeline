from __future__ import annotations

import numpy as np
import pytest

from jwst_stack.io import mask_scaled_sci
from jwst_stack.stack import sigma_clip_stack
from tests.conftest import gaussian_frame, make_wcs


def _frames(n: int = 3, spike_at=(16, 24)):
    wcs = make_wcs()
    base = gaussian_frame(wcs, peak_px=(16.0, 16.0), amp=100.0)
    frames = [base.copy() for _ in range(n)]
    frames[-1][spike_at] += 500.0
    return frames


def test_sigma_clip_removes_outlier():
    frames = _frames()
    stacked, coverage = sigma_clip_stack(frames, sigma=2.5, iterations=2)
    assert coverage[16, 24] == 3
    # Two clean frames (100 value is irrelevant here; the spike pixel base) are
    # both 0.0 (+gaussian wings ~0), so the clipped median must be ~0, not 500.
    assert stacked[16, 24] == pytest.approx(0.0, abs=5.0)


def test_coverage_counts_nan_frames():
    frames = _frames()
    frames[1][16, 16] = np.nan
    stacked, coverage = sigma_clip_stack(frames, sigma=3.0, iterations=2)
    assert coverage[16, 16] == 2
    # With 2 valid frames (< min for clipping) the plain median is used.
    assert np.isfinite(stacked[16, 16])


def test_fallback_small_overlap_two_frames():
    frames = _frames(n=2)
    stacked, coverage = sigma_clip_stack(frames, sigma=3.0, iterations=3)
    assert coverage.max() == 2
    assert np.isfinite(stacked[16, 16])


def test_mean_and_median_consistent_on_clean_pixels():
    frames = _frames(n=4)
    med, _ = sigma_clip_stack(frames, sigma=3.0, iterations=2, combine="median")
    mean, _ = sigma_clip_stack(frames, sigma=3.0, iterations=2, combine="mean")
    y, x = np.mgrid[0:64, 0:64]
    clean = (x - 10) ** 2 + (y - 10) ** 2 > 100
    assert np.allclose(med[clean], mean[clean], atol=1e-4, equal_nan=True)


def test_mask_scaled_sci_flags_dq_bit0():
    sci = np.ones((8, 8), dtype=np.float32)
    dq = np.zeros((8, 8), dtype=np.int64)
    dq[2, 2] = 0b1
    dq[3, 3] = 0b1100
    dq[4, 4] = 0b11
    out = mask_scaled_sci(sci, dq)
    assert np.isnan(out[2, 2])
    assert np.isnan(out[4, 4])
    assert out[3, 3] == 1.0
    assert np.isfinite(out).sum() == 64 - 2


def test_mask_scaled_sci_handles_nan_input():
    sci = np.ones((4, 4), dtype=np.float32)
    sci[1, 1] = np.nan
    sci[0, 0] = np.inf
    out = mask_scaled_sci(sci, dq=None)
    assert np.isnan(out[1, 1]) and np.isnan(out[0, 0])
    assert np.isfinite(out).sum() == 14