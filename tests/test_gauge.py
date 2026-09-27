"""Tests for the ``gauge`` step-2 subcommand.

The synthetic frames here are built so that the *measured* offset is known
exactly: the i2d and the frame share one WCS, and the frame's stars are drawn
at ``i2d_pixel + delta``.  Because the gauge converts frame pixels to sky and
then to i2d pixels, that is arithmetically identical to a header WCS that is
wrong by ``delta`` -- which is the thing being measured -- without needing a
reprojection.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from jwst_stack import gauge, io
from jwst_stack.gauge import GaugeError, GaugeGroup

from .conftest import make_wcs

LONGWAVE_CD = np.array([[0.063 / 3600.0, 0.0], [0.0, 0.063 / 3600.0]])
SHAPE = (64, 64)


def _star_pixels(n: int = 30) -> np.ndarray:
    """A deterministic, well-separated lattice of star positions."""
    side = int(np.ceil(np.sqrt(n)))
    xs = np.linspace(14, 50, side)
    grid = np.array([(x, y) for y in xs for x in xs])
    return grid[:n]


def _render(positions: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    ny, nx = SHAPE
    img = np.zeros(SHAPE, dtype=np.float32)
    yy, xx = np.mgrid[0:ny, 0:nx]
    sigma = 1.06  # ~2.5 px FWHM, the module's default
    for px, py in positions:
        img += 100.0 * np.exp(
            -((xx - px) ** 2 + (yy - py) ** 2) / (2.0 * sigma**2)
        )
    return img + rng.normal(0.0, 0.01, SHAPE).astype(np.float32)


def _write_i2d(path: Path, wcs, star_xy: np.ndarray) -> Path:
    rng = np.random.default_rng(7)
    sci = _render(star_xy, rng)
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(data=sci, header=wcs.to_header(), name="SCI"),
            fits.ImageHDU(data=np.ones(SHAPE, dtype=np.float32), name="WHT"),
        ]
    ).writeto(path, overwrite=True)
    return path


def _exposure(
    name: str,
    wcs,
    star_xy: np.ndarray,
    *,
    visit: int,
    detector: str,
    exposure: int = 1,
    pupil: str = "CLEAR",
    filter_name: str = "F335M",
    blank: bool = False,
) -> io.CalExposure:
    rng = np.random.default_rng(1000 + visit * 10 + len(detector))
    if blank:
        sci = rng.normal(0.0, 0.01, SHAPE).astype(np.float32)
    else:
        sci = _render(star_xy, rng)
    return io.CalExposure(
        path=Path(f"{name}.fits"),
        visit=visit,
        exposure=exposure,
        filter_name=filter_name,
        detector=detector,
        effexptm=1.0,
        wcs=wcs,
        shape=SHAPE,
        pupil=pupil,
        _sci=sci,
        _dq=np.zeros(SHAPE, dtype=np.int64),
    )


def _gauge_kwargs() -> dict:
    return dict(min_matched=8, max_i2d_stars=500, max_frame_stars=100, fwhm_px=2.5)


def test_measure_gauge_recovers_known_per_visit_offsets(tmp_path):
    """A known (+3, -2) px visit-2 offset must come back out of the gauge."""
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)

    exposures = [
        _exposure("v1a", wcs, stars, visit=1, detector="nrcalong"),
        _exposure("v1b", wcs, stars, visit=1, detector="nrcblong"),
        _exposure("v2a", wcs, stars + np.array([3.0, -2.0]), visit=2, detector="nrcalong"),
        _exposure("v2b", wcs, stars + np.array([3.0, -2.0]), visit=2, detector="nrcblong"),
    ]
    result = gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())

    assert set(result.groups) == {
        "v1_nrcalong", "v1_nrcblong", "v2_nrcalong", "v2_nrcblong"
    }
    v1 = result.visits[1]
    v2 = result.visits[2]
    assert v1.dx_median == pytest.approx(0.0, abs=0.15)
    assert v1.dy_median == pytest.approx(0.0, abs=0.15)
    assert v2.dx_median == pytest.approx(3.0, abs=0.15)
    assert v2.dy_median == pytest.approx(-2.0, abs=0.15)
    assert result.separation_from_anchor(2) == pytest.approx(np.hypot(3.0, 2.0), abs=0.2)


def test_anchor_is_the_visit_nearest_zero_not_the_lowest_number(tmp_path):
    """The whole point of measuring: the anchor need not be visit 1.

    This is the case where the ``min-visit`` heuristic would silently pin the
    wrong frame, so it is the case the subcommand exists to get right.  The
    offset stays inside the match radius (0.3 arcsec = 4.76 px at 0.063), as
    the real inter-visit pointing does; see
    ``test_offset_wider_than_the_match_radius_is_flagged_not_trusted`` for what
    happens outside it.
    """
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)

    exposures = [
        _exposure("v1", wcs, stars + np.array([2.0, 1.5]), visit=1, detector="nrcalong"),
        _exposure("v2", wcs, stars, visit=2, detector="nrcalong"),
    ]
    result = gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())

    assert result.anchor_visit == 2
    assert result.offset_median == pytest.approx(0.0, abs=0.2)
    assert result.visits[1].offset_median == pytest.approx(np.hypot(2.0, 1.5), abs=0.25)
    assert result.separation_from_anchor(1) == pytest.approx(np.hypot(2.0, 1.5), abs=0.25)


def test_offset_wider_than_the_match_radius_is_flagged_not_trusted(tmp_path):
    """A 7.8 px offset cannot be measured at a 4.76 px radius.

    Documenting a real limitation rather than a solved one.  The naive result
    is the dangerous kind: neighbouring stars still pair up, ``n_matched``
    clears the floor, and the returned median is 2.52 px with a *tight* scatter
    of 0.024 px -- indistinguishable, from the median and the scatter alone,
    from a genuine 2.5 px offset.  A lattice makes the wrong matches
    self-consistent, so no scatter-based test can catch this.

    What does catch it is that the median has grown to over half the search
    radius, which is where a real translation and a lucky pairing stop being
    separable.  The report says so instead of presenting the number as fact.
    """
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)

    exposures = [
        _exposure("v1", wcs, stars + np.array([6.0, 5.0]), visit=1, detector="nrcalong"),
        _exposure("v2", wcs, stars, visit=2, detector="nrcalong"),
    ]
    result = gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())

    # the misleading numbers, spelled out so a future change cannot hide them
    assert result.groups["v1_nrcalong"].n_matched > 8
    assert result.groups["v1_nrcalong"].offset_median == pytest.approx(2.5, abs=0.3)
    assert result.groups["v1_nrcalong"].offset_mad < 0.1  # deceptively tight

    assert gauge.suspect_groups(result) == ["v1_nrcalong"]
    text = gauge.format_gauge_report(result)
    assert "WARNING: 1 group(s) sit within" in text
    assert "not" in text and "distinguishable" in text
    assert "v1_nrcalong" in text


def test_match_rate_separates_a_real_offset_from_a_bad_match(tmp_path):
    """The diagnostic that actually works, where scatter does not.

    A real offset inside the radius matches nearly every star; a frame whose
    stars only partly overlap the i2d matches a fraction of them.
    ``offset_mad`` is tight in both cases, so the rate is the number to read -
    here it says the offset rests on 15 stars, not 30.
    """
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)

    # half the frame's stars have no counterpart in the i2d: displaced by more
    # than the 4.76 px match radius, so they are detected but never matched
    unmatched = stars + np.array([14.0, 14.0])
    mixed = np.vstack([stars, unmatched])

    exposures = [
        _exposure("v1", wcs, stars, visit=1, detector="nrcalong"),
        _exposure("v2", wcs, mixed, visit=2, detector="nrcalong"),
    ]
    result = gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())

    group = result.groups["v2_nrcalong"]
    # all 30 real counterparts matched, but the frame also offered stars with no
    # counterpart, so the rate is well below 1 and says how many backed it
    assert group.n_matched == 30
    assert group.n_frame_stars > group.n_matched
    assert 0.5 < group.match_rate < 1.0
    # the surviving matches still give the right answer
    assert group.dx_median == pytest.approx(0.0, abs=0.2)
    assert group.dy_median == pytest.approx(0.0, abs=0.2)
    # ...and the rate is on the record and in the report
    assert "rate" in gauge.format_gauge_report(result)
    assert "n_frame_stars" in result.as_groups_dict()["v2_nrcalong"]


def test_match_rate_is_nan_when_no_frame_stars_were_offered():
    assert GaugeGroup(1, "nrca1", 10, 0.0, 0.0, 0.0, 0.0).match_rate != (
        GaugeGroup(1, "nrca1", 10, 0.0, 0.0, 0.0, 0.0).match_rate
    )


def test_realistic_offsets_are_never_flagged(tmp_path):
    """The guard must not fire on the regime the shipped filters live in."""
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)

    exposures = [
        _exposure("v1", wcs, stars, visit=1, detector="nrcalong"),
        _exposure("v2", wcs, stars + np.array([0.14, -0.27]), visit=2,
                  detector="nrcalong"),
    ]
    result = gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())

    assert gauge.suspect_groups(result) == []
    assert "WARNING" not in gauge.format_gauge_report(result)


def test_offset_median_arcsec_uses_the_measured_i2d_scale(tmp_path):
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)
    exposures = [
        _exposure("v1", wcs, stars, visit=1, detector="nrcalong"),
        _exposure("v2", wcs, stars + np.array([4.0, 0.0]), visit=2, detector="nrcalong"),
    ]
    result = gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())

    assert result.i2d_scale_arcsec == pytest.approx(0.063, rel=1e-3)
    assert result.offset_median_arcsec == pytest.approx(
        result.offset_median * 0.063, rel=1e-6
    )


def test_measure_gauge_ignores_exposures_other_than_the_reference(tmp_path):
    """Only ``--exposure`` supplies stars, so other exposures must not add groups."""
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)

    exposures = [
        _exposure("e1", wcs, stars, visit=1, detector="nrcalong", exposure=1),
        _exposure("e2", wcs, stars + np.array([9.0, 9.0]), visit=1,
                  detector="nrcalong", exposure=2),
    ]
    result = gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())

    assert set(result.groups) == {"v1_nrcalong"}
    assert result.visits[1].offset_median == pytest.approx(0.0, abs=0.2)


def test_skipped_group_is_reported_and_excluded(tmp_path):
    """A group with no stars is listed, not silently averaged into the visit."""
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)

    exposures = [
        _exposure("v1a", wcs, stars, visit=1, detector="nrcalong"),
        _exposure("v1b", wcs, stars, visit=1, detector="nrcblong", blank=True),
    ]
    result = gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())

    assert set(result.groups) == {"v1_nrcalong"}
    assert result.visits[1].n_det == 1
    assert any("v1_nrcblong" in s for s in result.skipped)


def test_raises_when_no_group_clears_the_floor(tmp_path):
    """A gauge that cannot be read must not look like a gauge of zero."""
    wcs = make_wcs(cd=LONGWAVE_CD)
    blank = np.zeros(SHAPE, dtype=np.float32)
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, np.empty((0, 2)))

    exposures = [
        _exposure("v1", wcs, blank, visit=1, detector="nrcalong", blank=True)
    ]
    with pytest.raises(GaugeError):
        gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())


def test_raises_when_the_selection_matches_no_exposure_number(tmp_path):
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)
    exposures = [
        _exposure("e3", wcs, stars, visit=1, detector="nrcalong", exposure=3)
    ]
    with pytest.raises(GaugeError, match="exposure 1"):
        gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())


def test_summarise_visits_takes_the_median_not_the_mean():
    """One disagreeing detector must not drag the visit median."""
    groups = {
        "v1_nrca1": GaugeGroup(1, "nrca1", 500, 0.0, 0.0, 0.10, 0.01),
        "v1_nrca2": GaugeGroup(1, "nrca2", 5, 0.0, 0.0, 0.12, 0.01),
        "v1_nrcb1": GaugeGroup(1, "nrcb1", 500, 1.0, 0.0, 1.10, 0.01),
    }
    v = gauge.summarise_visits(groups)[1]

    assert v.dx_median == pytest.approx(0.0)  # median of (0, 0, 1)
    assert v.n_det == 3
    assert v.spread == pytest.approx(1.00)  # max - min of (0.10, 0.12, 1.10)
    assert v.internal_rms == pytest.approx(np.sqrt((0.0**2 + 0.0**2 + 1.0**2) / 3))


def test_summarise_visits_internal_rms_is_zero_for_identical_detectors():
    groups = {
        f"v1_det{i}": GaugeGroup(1, f"det{i}", 100, 0.4, -0.2, 0.45, 0.01)
        for i in range(4)
    }
    v = gauge.summarise_visits(groups)[1]
    assert v.internal_rms == pytest.approx(0.0)
    assert v.spread == pytest.approx(0.0)
    assert v.offset_median == pytest.approx(np.hypot(0.4, 0.2))


def test_format_gauge_report_states_the_anchor_and_the_next_command(tmp_path):
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)
    # the measured F335M v4 offset, so this exercises the realistic regime
    exposures = [
        _exposure("v1", wcs, stars, visit=1, detector="nrcalong"),
        _exposure("v2", wcs, stars + np.array([0.31, -0.38]), visit=2,
                  detector="nrcalong"),
    ]
    text = gauge.format_gauge_report(
        gauge.measure_gauge(exposures, i2d, **_gauge_kwargs()),
        filter_label="F335M",
    )

    assert "step 2 gauge probe filter=F335M" in text
    assert "measured scale 0.063" in text
    assert "anchor by per-visit median: visit 1" in text
    assert "raw inter-visit pointing error" in text
    assert "per (visit x detector) raw frame WCS vs official i2d" in text
    assert "nrcalong" in text
    assert "WARNING" not in text


def test_write_gauge_outputs_uses_the_flat_peer_schema(tmp_path):
    """The JSON must stay a peer of the existing *_step2_detectors.json files."""
    wcs = make_wcs(cd=LONGWAVE_CD)
    stars = _star_pixels()
    i2d = _write_i2d(tmp_path / "i2d.fits", wcs, stars)
    exposures = [
        _exposure("v1", wcs, stars, visit=1, detector="nrcalong"),
        _exposure("v2", wcs, stars, visit=2, detector="nrcalong"),
    ]
    result = gauge.measure_gauge(exposures, i2d, **_gauge_kwargs())

    json_out = tmp_path / "f335m_step2_detectors.json"
    log_out = tmp_path / "f335m_step2.log"
    gauge.write_gauge_outputs(
        result, json_path=json_out, log_path=log_out, filter_label="F335M"
    )

    import json

    payload = json.loads(json_out.read_text())
    assert set(payload) == {"v1_nrcalong", "v2_nrcalong"}
    # F187N's committed JSON has exactly visit/detector/n_matched/dx/dy/offset/
    # mad; this is that set plus the match-rate denominator, so it stays a peer.
    assert set(payload["v1_nrcalong"]) == {
        "visit", "detector", "n_matched", "dx_median", "dy_median",
        "offset_median", "offset_mad", "n_frame_stars",
    }
    assert payload["v1_nrcalong"]["n_frame_stars"] == 30
    assert "anchor by per-visit median" in log_out.read_text()


def test_filter_label_marks_a_narrowband_pupil():
    def make(pupil):
        return [
            io.CalExposure(Path("a"), 1, 1, "F444W", "nrcalong", 1.0, make_wcs(),
                           SHAPE, pupil)
        ]

    assert gauge.filter_label(make("CLEAR")) == "F444W"
    assert gauge.filter_label(make("F470N")) == "F444W;F470N"
    assert gauge.filter_label([]) == ""
