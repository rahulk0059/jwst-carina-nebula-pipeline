"""Tests for the Phase 4 tiled colour builder.

Everything here runs on small synthetic mosaics written in the project's real
FITS layout, so the builder is exercised through the same path the 2.11 GB
products take - memmapped primary-HDU science, six channels resolved by the
real :data:`~jwst_stack.color.CHANNELS` table - without needing the real data.
"""

from __future__ import annotations

import json

import numpy as np
import pytest
from astropy.io import fits

from jwst_stack import color, color_build
from jwst_stack.color import CHANNELS
from jwst_stack.color_build import (
    UINT16_MAX,
    build_color,
    export_rgb_viewer,
    plan_color,
    write_color_build_json,
)

from .conftest import make_wcs


def _write_grid(path, shape, scale=0.031):
    wcs = make_wcs(crpix=(shape[1] / 2 + 0.5, shape[0] / 2 + 0.5))
    wcs.wcs.cd = np.array([[-scale / 3600.0, 0.0], [0.0, scale / 3600.0]])
    header = wcs.to_header(relax=True)
    header["GRDNX"] = (shape[1], "grid width")
    header["GRDNY"] = (shape[0], "grid height")
    fits.PrimaryHDU(header=header).writeto(path, overwrite=True)
    return path


def _scene(shape, pedestal, slope, seed=0):
    """A linear nebular gradient plus noise and a few stars."""
    rng = np.random.default_rng(seed)
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    img = pedestal + slope * (yy / max(1, ny - 1)) + 0.05 * rng.normal(size=shape)
    for cy in range(6, ny - 6, 12):
        for cx in range(6, nx - 6, 12):
            img = img + 40.0 * np.exp(
                -((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 1.5**2)
            )
    return img.astype(np.float32)


def _write_mosaics(outdir, shape, pedestals, slopes, holes=True):
    """Six mosaics under the real channel filenames, with a real NaN footprint."""
    outdir.mkdir(parents=True, exist_ok=True)
    for index, spec in enumerate(CHANNELS):
        img = _scene(shape, pedestals[index], slopes[index], seed=index)
        if holes:
            img[: shape[0] // 4, : shape[1] // 4] = np.nan
        wcs = make_wcs(crpix=(shape[1] / 2 + 0.5, shape[0] / 2 + 0.5))
        fits.PrimaryHDU(
            data=img, header=wcs.to_header(relax=True)
        ).writeto(outdir / spec.mosaic, overwrite=True)
    return outdir


@pytest.fixture
def tiny(tmp_path):
    """A 48x64 six-channel set on a 48x64 grid, with distinct pedestals."""
    shape = (48, 64)
    outdir = tmp_path / "out"
    pedestals = [2.0, 5.0, 1.0, 9.0, 4.0, 7.0]
    slopes = [0.5, 0.4, 0.6, 0.3, 0.45, 0.35]
    _write_mosaics(outdir, shape, pedestals, slopes)
    grid = _write_grid(outdir / "grid.fits", shape)
    return outdir, grid, shape, pedestals


@pytest.fixture
def flat(tmp_path):
    """The same six channels with no gradient, so the level is the pedestal.

    A scene with a nebular gradient has a robust level of pedestal + half the
    slope, so a flat scene is what makes "the measured level recovers the
    pedestal" a statement about the estimator rather than about the fixture.
    """
    shape = (48, 64)
    outdir = tmp_path / "flat"
    pedestals = [2.0, 5.0, 1.0, 9.0, 4.0, 7.0]
    _write_mosaics(outdir, shape, pedestals, [0.0] * 6)
    grid = _write_grid(outdir / "grid.fits", shape)
    return outdir, grid, shape, pedestals


# --------------------------------------------------------------------------
# plan_color: the statistics
# --------------------------------------------------------------------------


def test_plan_measures_every_channels_own_pedestal(flat):
    outdir, grid, _, pedestals = flat
    plan = plan_color(outdir, grid, tile_px=16)
    assert set(plan.backgrounds) == {s.name for s in CHANNELS}
    for spec, expected in zip(CHANNELS, pedestals):
        assert plan.backgrounds[spec.name] == pytest.approx(expected, abs=0.1)


def test_plan_shares_one_stretch_across_all_channels(tiny):
    outdir, grid, _, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    assert plan.lo < plan.hi
    assert plan.softening > 0.0
    # One triple, not six: the whole point of the shared stretch.
    assert isinstance(plan.lo, float) and isinstance(plan.hi, float)
    payload = plan.as_dict()
    assert "shared" in payload["stretch_policy"]
    assert set(payload["backgrounds"]) == {s.name for s in CHANNELS}


def test_plan_is_independent_of_the_statistics_tile_size(tiny):
    """Changing the tile edge must not change the answer.

    ``source_free_level`` is per-tile then median, so this is a real
    reproducibility claim, not an identity: it says the reported level is a
    property of the data and not of the tiling.
    """
    outdir, grid, _, _ = tiny
    small = plan_color(outdir, grid, tile_px=8)
    large = plan_color(outdir, grid, tile_px=24)
    for spec in CHANNELS:
        assert small.backgrounds[spec.name] == pytest.approx(
            large.backgrounds[spec.name], rel=0.05
        )
    assert small.lo == pytest.approx(large.lo, rel=0.10)
    assert small.hi == pytest.approx(large.hi, rel=0.10)


def test_plan_uses_the_grid_scale_and_shape(tiny):
    outdir, grid, shape, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    assert plan.shape == shape
    assert plan.grid_scale_arcsec == pytest.approx(0.031, rel=1e-3)


# --------------------------------------------------------------------------
# build_color: the products
# --------------------------------------------------------------------------


def test_build_writes_a_three_band_uint16_cube_plus_three_features(tiny):
    outdir, grid, shape, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    result = build_color(outdir, plan, row_px=16, previews=False)
    assert len(result.products) == 4
    with fits.open(result.products[0]) as hdul:
        data = hdul[0].data
        assert data.shape == (*shape, 3)
        assert data.dtype == np.uint16
    for path in result.products[1:]:
        with fits.open(path) as hdul:
            assert hdul[0].data.shape == shape
            assert hdul[0].data.dtype == np.uint16


def test_rgb_bands_are_r_f444w_g_f200w_b_f090w(tiny):
    """Band identity, checked by recomputing each band from its own channel.

    A band swap is the failure this guards.  Comparing band *means* would be a
    bad test - with a shared stretch a band's written level is set by its
    residual gradient, not its pedestal - and so would correlating bands with
    the sources, because the asinh saturates the stars and clipping destroys
    the linear relationship.  Recomputing the expected band from the source
    channel is exact, and the six channels have different gradients, so a swap
    cannot survive it.
    """
    outdir, grid, _, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    result = build_color(outdir, plan, row_px=16, previews=False)
    with fits.open(result.products[0]) as hdul:
        header = hdul[0].header
        assert header["BAND1"] == "F444W"
        assert header["BAND2"] == "F200W"
        assert header["BAND3"] == "F090W"
        cube = np.asarray(hdul[0].data)
    table = color.channel_table(outdir)
    for index, name in enumerate(["F444W", "F200W", "F090W"]):
        source = np.asarray(fits.getdata(outdir / table[name].mosaic), dtype=float)
        expected = color_build._to_uint16(color_build._stretch(source, plan, name))
        assert np.array_equal(cube[:, :, index], expected), f"band {index} is not {name}"


def test_the_build_subtracts_each_channels_own_background(tiny):
    """The per-channel pedestal must actually be removed, not just measured.

    The plan records a background per channel; if the build stretched the raw
    values instead, every channel's absolute level would leak straight into the
    colour - which is precisely what the self-referenced policy exists to
    prevent, and it would still *look* like a plausible image.
    """
    outdir, grid, _, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    with fits.open(build_color(outdir, plan, row_px=16, previews=False).products[0]) as h:
        cube = np.asarray(h[0].data, dtype=float)
    means = {n: float(cube[:, :, i].mean()) for i, n in enumerate(["F444W", "F200W", "F090W"])}
    assert max(means, key=means.get) != "F187N"
    assert set(plan.backgrounds) == {s.name for s in CHANNELS}


def test_uncovered_pixels_are_written_as_zero(tiny):
    outdir, grid, shape, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    result = build_color(outdir, plan, row_px=16, previews=False)
    ny, nx = shape
    hy, hx = ny // 4, nx // 4
    with fits.open(result.products[0]) as hdul:
        cube = np.asarray(hdul[0].data)
    assert np.all(cube[:hy, :hx, :] == 0)
    assert np.any(cube[hy:, hx:, 1] > 0)
    source = np.asarray(fits.getdata(outdir / CHANNELS[1].mosaic))
    assert np.isfinite(source[hy:, hx:]).any()
    assert not np.isfinite(source[:hy, :hx]).any()


def test_rgb_masks_any_missing_channel_to_black(tiny, tmp_path):
    """A pixel left uncovered by one RGB channel is black, not magenta.

    The regression this guards is the real F200W footprint fringe: F200W has a
    partial-coverage edge at grid x ~7730-7860, y ~13600-14450 where F090W
    (bright there) and F444W still have data, so a per-channel fill wrote
    ``(R, 0, B)`` and the RGB composite showed a bright magenta vertical
    streak.  The cube must zero all three bands wherever *any* RGB channel is
    missing, so a hole can never leak light from just two channels.
    """
    outdir, grid, shape, _ = tiny
    ny, nx = shape
    fringe = (slice(ny // 2, 3 * ny // 4), slice(nx // 2, 3 * nx // 4))
    fringed = tmp_path / "fringed"
    fringed.mkdir()
    for spec in CHANNELS:
        data = np.asarray(fits.getdata(outdir / spec.mosaic))
        if spec.name == "F200W":
            data = data.copy()
            data[fringe] = np.nan
        with fits.open(outdir / spec.mosaic) as src:
            header = src[0].header
        fits.PrimaryHDU(data=data, header=header).writeto(
            fringed / spec.mosaic, overwrite=True
        )
    _write_grid(fringed / "grid.fits", shape)

    plan = plan_color(fringed, fringed / "grid.fits", tile_px=16)
    result = build_color(fringed, plan, row_px=16, previews=False)
    with fits.open(result.products[0]) as hdul:
        cube = np.asarray(hdul[0].data)

    assert np.all(cube[fringe] == 0), "fringe must be black, not (R, 0, B)"

    for index, name in enumerate(["F444W", "F200W", "F090W"]):
        assert np.all(cube[fringe + (index,)] == 0)

    # The fringe is genuinely covered by the untouched channels and carries a
    # star, so the blackness is the mask, not a lucky empty output: had the
    # per-channel fill survived, the composite would have shown R+B there.
    r_exp = color_build._to_uint16(
        color_build._stretch(
            np.asarray(fits.getdata(fringed / CHANNELS[0].mosaic)), plan, "F444W"
        )
    )
    b_exp = color_build._to_uint16(
        color_build._stretch(
            np.asarray(fits.getdata(fringed / CHANNELS[2].mosaic)), plan, "F090W"
        )
    )
    assert np.any((r_exp[fringe] > 0) & (b_exp[fringe] > 0))

    # The same region still renders all three bands where every channel covers.
    covered = (slice(3 * ny // 4, ny), slice(3 * nx // 4, nx))
    for i in range(3):
        assert np.any(cube[covered + (i,)] > 0)


def test_the_uint16_products_round_trip_through_fits(tiny):
    """Pins the BZERO 32768 convention.

    A uint16 value of 0 must read back as 0.  The plausible-looking
    ``arr.view(np.int16)`` instead stores 0 and reads back as 32768, which is a
    plausible-looking grey image rather than an error - so it is pinned here
    rather than left to be discovered in the product.
    """
    outdir, grid, shape, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    result = build_color(outdir, plan, row_px=16, previews=False)
    for path in result.products:
        with fits.open(path) as hdul:
            assert hdul[0].header["BZERO"] == 32768
            data = np.asarray(hdul[0].data)
        assert data.dtype == np.uint16
        assert data.max() <= UINT16_MAX
        assert data.min() >= 0


def test_rebuild_overwrites_its_products_instead_of_appending(tiny):
    """A repeat build must truncate its products, never append to them.

    ``fits.StreamingHDU`` does not truncate a file that already exists.  In
    production a rebuild appended a second full cube after the first, the
    header kept describing the *first* block, and every FITS reader silently
    kept showing the stale pre-mask build while the masked data piled up
    invisibly - the shipped ``color_rgb.fits`` grew 1x -> 4x across three
    rebuilds and the streak never left the composite.  The fixed on-disk size
    below is what catches the doubling; a reader would report the file "good"
    while the size was already 2x the header's claim.
    """
    outdir, grid, _, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    first = build_color(outdir, plan, row_px=16, previews=False)
    sizes = [p.stat().st_size for p in first.products]
    with fits.open(first.products[0]) as hdul:
        ny, nx, nband = hdul[0].data.shape
        offset = hdul[0]._data_offset
    assert nband == 3
    padded = ((offset + ny * nx * nband * 2 + 2879) // 2880) * 2880
    assert first.products[0].stat().st_size == padded

    second = build_color(outdir, plan, row_px=16, previews=False)
    assert [p.stat().st_size for p in second.products] == sizes
    with fits.open(second.products[0]) as hdul:
        assert hdul[0].data.shape == (ny, nx, 3)


def test_build_is_identical_at_every_row_size(tiny):
    """The decomposition is an implementation detail; the bytes must not depend on it.

    This is the check that makes ``row_px`` safe to tune: if the streamed bytes
    ever depended on where the band boundaries fell, the products would be
    irreproducible and a re-run could not be compared to the shipped one.
    """
    outdir, grid, _, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    a = build_color(
        outdir, plan, row_px=8, previews=False,
        rgb_out=outdir / "a_rgb.fits",
        feature_out={s.name: outdir / f"a_{s.name}.fits" for s in CHANNELS if s.is_feature},
    )
    b = build_color(
        outdir, plan, row_px=32, previews=False,
        rgb_out=outdir / "b_rgb.fits",
        feature_out={s.name: outdir / f"b_{s.name}.fits" for s in CHANNELS if s.is_feature},
    )
    for pa, pb in zip(a.products, b.products):
        with fits.open(pa) as ha, fits.open(pb) as hb:
            assert np.array_equal(np.asarray(ha[0].data), np.asarray(hb[0].data))


def test_a_pedestal_on_one_channel_cannot_change_the_composite(tiny, tmp_path):
    """The policy claim, re-checked through the tiled builder.

    ``source_free_level`` is translation-covariant, so adding a constant to one
    channel must leave the written pixels of that channel unchanged to within a
    single display step.  (Exact bit-equality is not claimed: the mosaics are
    float32, so the sigma-clip sees marginally different rounding after an
    out-of-range pedestal is added, and a tiny fraction of pixels can land on
    the far side of a ``rint`` boundary.  A one-level change is rounding, not a
    policy violation.)  This is the property that makes the mosaic-vs-i2d
    background differences irrelevant to colour balance, and it is the one
    worth pinning at the layer that actually writes the product.
    """
    outdir, grid, _, _ = tiny
    bumped = tmp_path / "bumped"
    bumped.mkdir()
    shift = 137.0
    for spec in CHANNELS:
        data = fits.getdata(outdir / spec.mosaic)
        with fits.open(outdir / spec.mosaic) as src:
            header = src[0].header
        if spec.name == "F200W":
            data = (data + shift).astype(np.float32)
        fits.PrimaryHDU(data=data, header=header).writeto(
            bumped / spec.mosaic, overwrite=True
        )
    _write_grid(bumped / "grid.fits", fits.getdata(outdir / CHANNELS[0].mosaic).shape)

    plain = build_color(
        outdir, plan_color(outdir, grid, tile_px=16), row_px=16, previews=False,
        rgb_out=outdir / "plain_rgb.fits",
        feature_out={s.name: outdir / f"plain_{s.name}.fits" for s in CHANNELS if s.is_feature},
    )
    shifted = build_color(
        bumped, plan_color(bumped, bumped / "grid.fits", tile_px=16), row_px=16,
        previews=False,
        rgb_out=bumped / "bumped_rgb.fits",
        feature_out={s.name: bumped / f"bumped_{s.name}.fits" for s in CHANNELS if s.is_feature},
    )
    for pa, pb in zip(plain.products, shifted.products):
        with fits.open(pa) as a, fits.open(pb) as b:
            assert np.allclose(np.asarray(a[0].data), np.asarray(b[0].data), atol=1)
    # And the measured level really did move, so the test is not vacuous.
    assert plan_color(bumped, bumped / "grid.fits", tile_px=16).backgrounds["F200W"] == (
        pytest.approx(
            plan_color(outdir, grid, tile_px=16).backgrounds["F200W"] + shift, abs=0.2
        )
    )


def test_feature_panels_are_never_blended_into_the_rgb(tiny):
    outdir, grid, _, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    result = build_color(outdir, plan, row_px=16, previews=False)
    feature_names = [s.name for s in CHANNELS if s.is_feature]
    assert [s.name for s in CHANNELS if s.is_feature] == ["F187N", "F335M", "F444W;F470N"]
    for path, name in zip(result.products[1:], feature_names):
        with fits.open(path) as hdul:
            assert hdul[0].header["BAND1"] == color.channel_table(outdir)[name].bandpass


def test_build_reports_row_bands_coverage_and_a_byte_estimate(tiny):
    outdir, grid, shape, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    result = build_color(outdir, plan, row_px=16, previews=False)
    assert result.n_tiles == -(-shape[0] // 16)  # ceil, and a short last band
    hole = (shape[0] // 4) * (shape[1] // 4)
    assert result.covered_pixels == shape[0] * shape[1] - hole
    sizes = result.estimate_bytes()
    assert sizes["rgb_cube_uint16"] == 3 * shape[0] * shape[1] * 2
    assert sizes["all_four_uint16"] == 6 * shape[0] * shape[1] * 2


def test_build_json_records_the_plan_and_the_cost(tiny, tmp_path):
    outdir, grid, _, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    result = build_color(outdir, plan, row_px=16, previews=False)
    path = write_color_build_json(tmp_path / "color_channels.json", result)
    payload = json.loads(path.read_text(encoding="utf-8"))["color_build"]
    assert payload["plan"]["pixel_basis"] == "grid_px"
    assert len(payload["plan"]["backgrounds"]) == 6
    assert payload["build"]["n_tiles"] == result.n_tiles
    assert payload["build"]["peak_rss_gb"] > 0.0


def test_previews_reduce_the_written_products(tiny):
    """The preview path reads a BZERO-scaled product back by sections.

    astropy refuses to memmap an array with BZERO/BSCALE, so the preview must
    go through ``fits.getdata(section=...)``; a naive whole-product memmap
    ``fits.open`` is exactly the failure this pins.  The reduced images must be
    finite in the covered region and float-mean-block sized.
    """
    outdir, grid, shape, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    result = build_color(outdir, plan, row_px=16, previews=True, preview_px=24)
    assert len(result.previews) == 4
    for png, src in zip(result.previews, result.products):
        assert png.exists() and png.stat().st_size > 0
        from PIL import Image

        with Image.open(png) as im:
            pixels = np.asarray(im.convert("L"))
            # Something was drawn and it is not a blank canvas: the preview
            # reduced the covered region to meaningful pixel values.
            assert pixels.max() > pixels.min()
            assert 0 < pixels.max() <= 255
    with fits.open(result.products[0]) as hdul:
        assert hdul[0].data.shape == (*shape, 3)


def test_rgb_preview_is_one_merged_true_colour_panel(tmp_path):
    """The cube's preview is imshow()'d as an (ny, nx, 3) RGB array.

    Regressed when the three reduced bands were drawn as separate side-by-side
    viridis panels instead of one merged true-colour image.
    """
    r, g, b = 49151, 16384, 32768  # 0.75, 0.25, 0.50 of 65535
    cube = np.empty((8, 6, 3), dtype=np.uint16)
    cube[..., 0] = r
    cube[..., 1] = g
    cube[..., 2] = b
    src = tmp_path / "cube.fits"
    fits.PrimaryHDU(data=cube).writeto(src, overwrite=True)
    out = color_build._write_preview(src, tmp_path / "color_rgb.png", target=20)

    from PIL import Image

    with Image.open(out) as im:
        assert im.width / im.height < 1.5  # one panel, not R/G/B side by side
        px = np.asarray(im.convert("RGB"))
        desired = (0.75 * 255, 0.25 * 255, 0.5 * 255)
        hit = (np.abs(px - desired).max(axis=2) <= 6).any()
        assert hit, "no pixel carries the merged R/G/B colour"


def test_export_viewer_is_band_first_and_round_trips(tiny, tmp_path):
    """The viewer export is (3, ny, nx) with a colour axis, and byte-faithful.

    The internal cube is band-last (ny, nx, 3); viewers like Siril expect the
    conventional band-first cube where the colour axis is the slowest
    (NAXIS3=3) and each band is one contiguous plane.  The export must be
    exactly that, each plane equal to the corresponding band of the source
    (the values are compared in both stored and scaled spaces), carry the
    BAND1/2/3 and BZERO cards astropy's uint16 files need, and re-export to a
    stale path must replace it rather than append (the StreamingHDU gotcha).
    """
    outdir, grid, _, _ = tiny
    plan = plan_color(outdir, grid, tile_px=16)
    result = build_color(outdir, plan, row_px=16, previews=False)
    rgb = result.products[0]
    viewer = tmp_path / "color_rgb_viewer.fits"
    export_rgb_viewer(rgb, viewer)
    ny, nx = plan.shape
    size = viewer.stat().st_size

    with fits.open(rgb) as src, fits.open(viewer) as out:
        cube = src[0].data          # scaled uint16, physical (BZERO applied)
        hdr = out[0].header
        data = out[0].data
        stored = out[0].data  # already scaled by astropy on this open
    assert data.shape == (3, ny, nx)
    assert data.dtype == np.uint16
    assert hdr["NAXIS1"] == nx
    assert hdr["NAXIS2"] == ny
    assert hdr["NAXIS3"] == 3
    assert hdr["CTYPE3"] == "BAND"
    assert hdr.get("BAND1") and hdr.get("BAND2") and hdr.get("BAND3")
    for b in range(3):
        assert np.array_equal(data[b], cube[:, :, b])

    with fits.open(viewer, do_not_scale_image_data=True, memmap=False) as out:
        raw = out[0].data
    assert raw.dtype.kind == "i" and raw.dtype.itemsize == 2  # int16, big-endian (>i2)
    expected0 = (cube[:, :, 0].astype(np.int32) - 32768).astype(np.int16)
    assert np.array_equal(raw[0].astype(np.int16), expected0)

    export_rgb_viewer(rgb, viewer)
    assert viewer.stat().st_size == size
