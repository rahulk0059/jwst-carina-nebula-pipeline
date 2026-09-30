from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from jwst_stack.starcat import (
    brightest_separated,
    coverage_probe,
    covered_tiles,
    detect_tiled,
    open_mosaic,
    tile_coverage_fraction,
)

from .conftest import make_wcs


def _write_mosaic(path, data, coverage=None, wcs=None):
    """A mosaic in the project's layout: SCI in the primary HDU.

    *coverage* is written the way ``mosaic`` writes it - unsigned int16 depth
    with ``BZERO 32768`` - because that is the layout that breaks memmapping.
    """
    wcs = wcs or make_wcs(crpix=(data.shape[1] / 2 + 0.5, data.shape[0] / 2 + 0.5))
    hdus = [fits.PrimaryHDU(data=data.astype(np.float32), header=wcs.to_header(relax=True))]
    if coverage is not None:
        hdus.append(
            fits.ImageHDU(
                data=coverage.astype(np.uint16),
                name="COVERAGE",
                header=fits.Header({"BZERO": 32768, "BSCALE": 1}),
            )
        )
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


def _star_field(shape=(256, 256), sigma=1.5, amp=800.0, pedestal=10.0, seed=0):
    """A flat background with well-separated Gaussian stars."""
    rng = np.random.default_rng(seed)
    img = np.full(shape, pedestal, dtype=float) + rng.normal(0.0, 1.0, shape)
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    for cy in range(30, shape[0] - 30, 40):
        for cx in range(30, shape[1] - 30, 40):
            img += amp * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sigma**2))
    return img


# --------------------------------------------------------------------------
# open_mosaic: the two FITS-layout traps, in one place
# --------------------------------------------------------------------------


def test_open_mosaic_reads_science_from_the_primary_hdu(tmp_path):
    data = np.arange(64 * 64, dtype=np.float32).reshape(64, 64)
    path = _write_mosaic(tmp_path / "m.fits", data)
    opened = open_mosaic(path)
    assert opened.shape == (64, 64)
    assert np.array_equal(np.asarray(opened.sci), data)


def test_open_mosaic_survives_the_unsigned_coverage_extension(tmp_path):
    """COVERAGE is BITPIX 16 with BZERO 32768, so astropy refuses to memmap it.

    Reading it as the coverage source would raise on every real mosaic, so the
    probe is derived from the science array's NaN mask instead.
    """
    data = np.full((64, 64), np.nan, dtype=np.float32)
    data[16:48, 16:48] = 2.0
    path = _write_mosaic(tmp_path / "m.fits", data, coverage=np.zeros((64, 64), dtype=int))
    opened = open_mosaic(path)
    assert opened.coverage(16, 48, 16, 48) == 1.0
    assert opened.coverage(0, 16, 0, 16) == 0.0


def test_coverage_probe_reads_coverage_depth_greater_than_zero():
    cov = np.zeros((16, 16), dtype=np.int16)
    cov[:8] = 3
    probe = coverage_probe(coverage=cov)
    assert probe(0, 16, 0, 16) == pytest.approx(0.5)


def test_coverage_probe_from_science_counts_finite_pixels():
    sci = np.full((16, 16), np.nan)
    sci[:8] = 1.0
    assert coverage_probe(sci=sci)(0, 16, 0, 16) == pytest.approx(0.5)


def test_coverage_probe_prefers_the_explicit_coverage_array():
    sci = np.zeros((16, 16))
    cov = np.zeros((16, 16), dtype=np.int16)
    cov[:4, :4] = 1
    # Sampled on an 8 px stride, so only (0, 0) of the covered corner is seen.
    assert coverage_probe(sci=sci, coverage=cov)(0, 16, 0, 16) == pytest.approx(0.25)


def test_coverage_probe_of_nothing_is_none():
    assert coverage_probe() is None


def test_open_mosaic_rejects_a_file_with_no_primary_data(tmp_path):
    path = tmp_path / "empty.fits"
    fits.HDUList([fits.PrimaryHDU()]).writeto(path, overwrite=True)
    with pytest.raises(ValueError, match="no science data"):
        open_mosaic(path)


# --------------------------------------------------------------------------
# tile selection
# --------------------------------------------------------------------------


def test_tile_coverage_fraction_uses_depth_greater_than_zero():
    cov = np.zeros((16, 16), dtype=np.int16)
    cov[:8] = 2
    assert tile_coverage_fraction(cov, 0, 16, 0, 16) == pytest.approx(0.5)


def test_covered_tiles_drops_half_empty_tiles():
    cov = np.zeros((64, 64), dtype=np.int16)
    cov[:32, :] = 1
    tiles = covered_tiles((64, 64), coverage_probe(coverage=cov), 32)
    assert tiles == [(0, 32, 0, 32), (0, 32, 32, 64)]


def test_covered_tiles_without_a_probe_returns_every_tile():
    assert len(covered_tiles((64, 64), None, 32)) == 4


def test_covered_tiles_spreads_a_capped_sample_across_the_field():
    cov = np.ones((640, 640), dtype=np.int16)
    tiles = covered_tiles((640, 640), coverage_probe(coverage=cov), 160, n_tiles=3)
    assert len(tiles) == 3
    assert tiles == sorted(tiles), "the sample should stay in row-major order"
    assert tiles[-1][0] > tiles[0][0], "the sample should not sit in one row"


def test_covered_tiles_cap_above_the_available_count_returns_all():
    cov = np.ones((64, 64), dtype=np.int16)
    probe = coverage_probe(coverage=cov)
    assert len(covered_tiles((64, 64), probe, 32, n_tiles=99)) == 4


# --------------------------------------------------------------------------
# detect_tiled
# --------------------------------------------------------------------------


def test_detect_tiled_finds_the_stars_in_a_field():
    img = _star_field()
    x, y, _ = detect_tiled(img, covered_tiles(img.shape, None, 128), halo_px=8)
    assert x.size >= 10
    assert np.all(np.isfinite(x)) and np.all(np.isfinite(y))


def test_a_star_in_a_halo_wide_enough_to_swallow_two_tiles_is_counted_once():
    """The exact property the interior test exists for.

    A 32 px star at (40, 40) with 32 px halos is visible from three different
    tiles.  If the halo band were kept, it would be counted three times.
    """
    img = np.full((128, 128), 10.0)
    yy, xx = np.mgrid[0:128, 0:128]
    img += 900.0 * np.exp(-((xx - 40) ** 2 + (yy - 40) ** 2) / (2 * 1.5**2))
    for tile_px in (32, 64, 128):
        fx, fy, _ = detect_tiled(img, covered_tiles(img.shape, None, tile_px), halo_px=32)
        real = (np.abs(fx - 40) < 1.0) & (np.abs(fy - 40) < 1.0)
        assert real.sum() == 1, f"tiled at {tile_px} px it was counted {real.sum()} times"


def test_detect_tiled_does_not_miss_a_star_by_tiling_it_finely():
    """Every star of the coarse run is found again on the fine grid.

    Equality is not expected, only coverage: DAOStarFinder's threshold is the
    MAD of the patch it is handed, so re-tiling changes the noise estimate
    slightly and can add or drop a faint edge detection.  What must not happen
    is a star going missing, or one star being counted twice.
    """
    img = _star_field(shape=(256, 256))
    fx, fy, _ = detect_tiled(img, covered_tiles(img.shape, None, 64), halo_px=16)
    cx, cy, _ = detect_tiled(img, covered_tiles(img.shape, None, 256), halo_px=16)
    assert _match_count((fx, fy), (cx, cy)) == cx.size
    assert fx.size <= cx.size + 4, "the fine grid inflated the count: double counting"
    assert _match_count((cx, cy), (fx, fy)) >= 25, "the fine grid lost real stars"


def _match_count(found, expected, radius: float = 0.5) -> int:
    """How many of *found* sit within *radius* of some point in *expected*."""
    fx, fy = found
    ex, ey = expected
    if fx.size == 0:
        return 0
    d2 = (fx[:, None] - ex[None, :]) ** 2 + (fy[:, None] - ey[None, :]) ** 2
    return int(np.count_nonzero(np.min(d2, axis=1) <= radius**2))


def test_detect_tiled_keep_margin_keeps_stars_away_from_the_expanded_edge():
    img = _star_field(shape=(128, 128))
    tiles = [(0, 64, 0, 64)]
    x, _, _ = detect_tiled(img, tiles, halo_px=8, keep_margin=12.0)
    assert np.all(x >= 12) and np.all(x < 64 + 8 - 12)


def test_detect_tiled_max_stars_bounds_the_catalogue():
    img = _star_field(shape=(256, 256))
    x, _, _ = detect_tiled(img, covered_tiles(img.shape, None, 128), halo_px=8, max_stars=3)
    assert x.size == 3


def test_detect_tiled_returns_empty_on_a_blank_field():
    blank = np.zeros((64, 64), dtype=float)
    x, y, f = detect_tiled(blank, covered_tiles(blank.shape, None, 64))
    assert x.size == y.size == f.size == 0


def test_detect_tiled_reports_progress():
    img = _star_field(shape=(128, 128))
    seen = []
    detect_tiled(
        img,
        covered_tiles(img.shape, None, 64),
        halo_px=8,
        progress=lambda i, n: seen.append((i, n)),
    )
    assert seen and seen[-1][0] == seen[-1][1] == 4


# --------------------------------------------------------------------------
# brightest_separated
# --------------------------------------------------------------------------


def test_brightest_separated_keeps_the_bright_member_of_a_close_pair():
    x = np.array([10.0, 12.0, 50.0])
    y = np.array([10.0, 10.0, 50.0])
    flux = np.array([100.0, 900.0, 500.0])
    picks = brightest_separated(x, y, flux, 3, separation=5.0)
    assert set(picks.tolist()) == {1, 2}


def test_brightest_separated_respects_the_count():
    x = np.arange(10.0)
    y = np.zeros(10)
    flux = np.arange(10.0)
    assert brightest_separated(x, y, flux, 3, separation=2.0).size == 3


def test_brightest_separated_of_nothing_is_empty():
    empty = np.empty(0)
    assert brightest_separated(empty, empty, empty, 5, 1.0).size == 0
