"""Tests for pupil-aware selection.

The bug this exists to prevent: ``FILTER`` is ``F444W`` for *both* the CLEAR and
the F470N exposure of program 2731, and selection filtered on detector and
filter only.  ``--filter F444W`` therefore matched all 80 files, two bandpasses
and two different science goals, and reported success throughout.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from jwst_stack import cli, download, io

from .test_download import _make_row, _products_table

SHAPE = (8, 8)
PREFIX = "jw02731001001_02103_00001_nrcalong_cal"


def _wcs() -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ["deg", "deg"]
    w.wcs.crpix = [4.5, 4.5]
    w.wcs.crval = [160.0, -58.0]
    w.wcs.cd = np.array([[-8.6e-6, 0.0], [0.0, 8.6e-6]])
    return w


def _write_cal(path: Path, filter_name: str, pupil: str, visit: int = 1) -> Path:
    sci = fits.Header()
    sci.update(_wcs().to_header())
    fits.HDUList(
        [
            fits.PrimaryHDU(
                header=fits.Header({"FILTER": filter_name, "PUPIL": pupil,
                                    "VISIT": visit, "EXPOSURE": 1})
            ),
            fits.ImageHDU(data=np.zeros(SHAPE, dtype=np.float32), header=sci, name="SCI"),
            fits.ImageHDU(data=np.zeros(SHAPE, dtype=np.int32), name="DQ"),
        ]
    ).writeto(path, overwrite=True)
    return path


def _exposure(filter_name: str, pupil: str, visit: int = 1, detector: str = "nrcalong"):
    return io.CalExposure(
        path=Path(f"{filter_name}-{pupil}-v{visit}-{detector}.fits"),
        visit=visit,
        exposure=1,
        filter_name=filter_name,
        detector=detector,
        effexptm=1.0,
        wcs=_wcs(),
        shape=SHAPE,
        pupil=pupil,
    )


def _f444w_pair() -> list[io.CalExposure]:
    """One CLEAR and one F470N frame: same FILTER, different bandpass."""
    return [
        _exposure("F444W", "CLEAR", visit=1, detector="nrcalong"),
        _exposure("F444W", "F470N", visit=1, detector="nrcalong"),
    ]


# --- header reading ---------------------------------------------------------


def test_read_cal_exposure_reads_the_pupil(tmp_path):
    """PUPIL lives in the primary header, next to FILTER, not in the SCI."""
    path = _write_cal(tmp_path / "clear_cal.fits", "F444W", "CLEAR")
    assert io.read_cal_exposure(path).pupil == "CLEAR"

    path = _write_cal(tmp_path / "grism_cal.fits", "F444W", "F470N")
    assert io.read_cal_exposure(path).pupil == "F470N"


def test_pupil_defaults_to_clear_when_the_header_lacks_it(tmp_path):
    """Old files, and short-wave files, have no PUPIL; they must still load."""
    path = _write_cal(tmp_path / "nopupil_cal.fits", "F335M", "")
    assert io.read_cal_exposure(path).pupil == "CLEAR"


def test_scan_exposures_sees_the_pupil(tmp_path):
    # find_cal_files globs *_cal.fits, so the visit/exposure goes before _cal
    _write_cal(tmp_path / "jw02731001001_02103_00001_nrcalong_cal.fits", "F444W", "CLEAR")
    _write_cal(tmp_path / "jw02731001002_02103_00003_nrcalong_cal.fits", "F444W", "F470N")
    found = {e.pupil for e in io.iter_exposures(tmp_path)}
    assert found == {"CLEAR", "F470N"}


# --- selection helpers ------------------------------------------------------


def test_mixed_pupils_reports_all_of_them():
    assert io.mixed_pupils(_f444w_pair()) == ["CLEAR", "F470N"]
    assert io.mixed_pupils([_exposure("F335M", "CLEAR")]) == ["CLEAR"]


def test_select_exposures_narrows_to_one_pupil():
    kept = io.select_exposures(_f444w_pair(), filter_name=["F444W"], pupil=["CLEAR"])
    assert [e.pupil for e in kept] == ["CLEAR"]


def test_select_exposures_defaults_to_every_pupil():
    """Bare filters keep their old meaning, so existing plans do not change."""
    kept = io.select_exposures(_f444w_pair(), filter_name=["F444W"])
    assert {e.pupil for e in kept} == {"CLEAR", "F470N"}


def test_select_exposures_returns_nothing_for_an_unknown_pupil():
    """A typo must not fall through to "every pupil"; the caller reports it."""
    assert io.select_exposures(_f444w_pair(), filter_name=["F444W"], pupil=["F090W"]) == []


# --- the CLI guard ----------------------------------------------------------


def test_cli_rejects_a_mixed_pupil_selection():
    """The regression test for the reported bug.

    Before ``--pupil`` existed this returned both frames and the caller happily
    stacked two bandpasses.  It must now refuse and name the way out.
    """
    with pytest.raises(SystemExit) as excinfo:
        cli._select_exposures(_f444w_pair(), ["nrcalong"], ["F444W"], None)
    message = str(excinfo.value)
    assert "CLEAR" in message and "F470N" in message
    assert "--pupil" in message
    assert "mix bandpasses" in message


def test_cli_accepts_a_mixed_tree_when_the_pupil_is_named():
    kept = cli._select_exposures(_f444w_pair(), ["nrcalong"], ["F444W"], ["CLEAR"])
    assert [e.pupil for e in kept] == ["CLEAR"]

    kept = cli._select_exposures(_f444w_pair(), ["nrcalong"], ["F444W"], ["F470N"])
    assert [e.pupil for e in kept] == ["F470N"]


def test_cli_does_not_require_a_pupil_for_a_single_pupil_filter():
    """F335M and the three short-wave filters must keep working unpupiled."""
    for filter_name in ("F090W", "F187N", "F200W", "F335M"):
        exposures = [_exposure(filter_name, "CLEAR", visit=v) for v in (1, 2)]
        assert len(cli._select_exposures(exposures, ["nrcalong"], [filter_name], None)) == 2


def test_cli_rejects_a_multi_pupil_request_that_is_itself_mixed():
    """``--pupil CLEAR F470N`` asks for the same mix and must be refused."""
    with pytest.raises(SystemExit, match="mix bandpasses"):
        cli._select_exposures(_f444w_pair(), ["nrcalong"], ["F444W"], ["CLEAR", "F470N"])


def test_narrowband_gauge_writes_to_its_own_file(tmp_path, monkeypatch, capsys):
    """F444W;F470N must not overwrite the F444W:CLEAR gauge.

    Both are ``FILTER=F444W``, so a stem built from the filter alone would have
    the second run silently clobber the first.
    """
    def measure(exposures, i2d, **kwargs):
        from jwst_stack.gauge import GaugeGroup, GaugeResult, summarise_visits

        groups = {
            "v1_nrcalong": GaugeGroup(1, "nrcalong", 30, 0.0, 0.0, 0.0, 0.0, 30)
        }
        result = GaugeResult(groups=groups, i2d_scale_arcsec=0.063)
        result.visits = summarise_visits(groups)
        result.anchor_visit = 1
        return result

    monkeypatch.setattr(cli.gauge_mod, "measure_gauge", measure)
    monkeypatch.setattr(
        cli,
        "_load_selected",
        lambda args: [e for e in _f444w_pair() if e.pupil == args.pupil[0]],
    )

    for pupil in ("CLEAR", "F470N"):
        args = cli.build_parser().parse_args(
            ["gauge", "--i2d", "x.fits", "--outdir", str(tmp_path),
             "--pupil", pupil]
        )
        cli.run_gauge(args)
    capsys.readouterr()

    written = sorted(p.name for p in tmp_path.iterdir())
    assert "f444w_step2_detectors.json" in written
    assert "f444w_f470n_step2_detectors.json" in written


def test_every_selection_command_accepts_the_pupil_flag():
    """A flag on one subparser and not the others is the `--stage2-kinds` trap."""
    for command in ("inspect", "group", "stack", "mosaic", "gauge"):
        extra = ["--i2d", "i2d.fits"] if command == "gauge" else []
        args = cli.build_parser().parse_args(
            [command, *extra, "--filter", "F444W", "--pupil", "CLEAR"]
        )
        assert args.pupil == ["CLEAR"], command


# --- download planning ------------------------------------------------------


@pytest.mark.parametrize(
    "spec,expected",
    [
        ("F444W", ("F444W", None)),
        ("f444w", ("F444W", None)),
        ("F444W;CLEAR", ("F444W", "CLEAR")),
        ("F444W;F470N", ("F444W", "F470N")),
        ("  f335m ; clear ", ("F335M", "CLEAR")),
    ],
)
def test_parse_filter_spec(spec, expected):
    assert download.parse_filter_spec(spec) == expected


def _stage2_table():
    """F444W in both pupils, two visits each, plus one combined i2d per pupil.

    Reuses ``test_download``'s table helpers so the filenames follow the same
    convention as the real MAST listing; ``parse_product_name`` has to fail on
    the combined i2d names for ``build_i2d_items`` to work at all.
    """
    rows = [
        _make_row("jw02731001001_02105_00001_nrcalong_cal.fits", "CAL", "F444W;CLEAR", 100),
        _make_row("jw02731001002_02105_00003_nrcalong_cal.fits", "CAL", "F444W;CLEAR", 100),
        _make_row("jw02731001001_02105_00001_nrcblong_cal.fits", "CAL", "F444W;F470N", 100),
        _make_row("jw02731001002_02105_00003_nrcblong_cal.fits", "CAL", "F444W;F470N", 100),
        _make_row("jw02731-o001_t017_nircam_clear-f444w_i2d.fits", "I2D", "F444W;CLEAR", 200),
        _make_row("jw02731-o001_t017_nircam_f444w-f470n_i2d.fits", "I2D", "F444W;F470N", 200),
        # a third filter, to prove the pupil filter does not leak across filters
        _make_row("jw02731001003_02105_00001_nrcblong_cal.fits", "CAL", "F335M", 100),
    ]
    return _products_table(rows)


def test_plan_stage2_narrows_to_the_named_pupil(tmp_path):
    prods = _stage2_table()
    plan = download.plan_stage2(
        prods, tmp_path / "data", tmp_path / "i2d", ["F444W;CLEAR"]
    )

    assert plan["filters"] == ["F444W;CLEAR"]
    assert plan["pupils"] == ["CLEAR"]
    assert len(plan["all_cal"]) == 2
    # every planned cal is CLEAR, and the combined i2d is the CLEAR one
    text = download.format_stage2_plan(plan)
    assert "F470N" not in text


def test_plan_stage2_bare_filter_still_means_every_pupil(tmp_path):
    """A bare name must not silently narrow to CLEAR and halve the plan."""
    prods = _stage2_table()
    plan = download.plan_stage2(prods, tmp_path / "data", tmp_path / "i2d", ["F444W"])

    assert plan["pupils"] == []
    assert len(plan["all_cal"]) == 4  # both pupils
    assert len(plan["combined"]) == 2


def test_bare_plan_warns_that_it_mixes_pupils(tmp_path):
    """A bare ``F444W`` must not quietly plan two bandpasses.

    The usual way to land here is a shell eating the ``;`` in ``F444W;CLEAR``,
    so the warning quotes the correct command.
    """
    prods = _stage2_table()
    text = download.format_stage2_plan(
        download.plan_stage2(prods, tmp_path / "data", tmp_path / "i2d", ["F444W"])
    )
    assert "WARNING: this plan mixes 2 pupils (CLEAR, F470N)" in text
    assert "QUOTE IT" in text
    assert 'F444W;CLEAR' in text


def test_narrow_plan_does_not_warn(tmp_path):
    prods = _stage2_table()
    text = download.format_stage2_plan(
        download.plan_stage2(prods, tmp_path / "data", tmp_path / "i2d",
                             ["F444W;CLEAR"])
    )
    assert "WARNING" not in text


def test_plan_stage2_keeps_itself_consistent_across_both_pupil_specs(tmp_path):
    """The two narrow plans must partition the bare plan, not overlap or miss."""
    prods = _stage2_table()
    args = (prods, tmp_path / "data", tmp_path / "i2d")
    bare = download.plan_stage2(*args, ["F444W"])
    clear = download.plan_stage2(*args, ["F444W;CLEAR"])
    grism = download.plan_stage2(*args, ["F444W;F470N"])

    def names(plan, key):
        return sorted(getattr(i, "filename") for i in plan[key])

    assert sorted(names(clear, "all_cal") + names(grism, "all_cal")) == names(bare, "all_cal")
    assert sorted(names(clear, "combined") + names(grism, "combined")) == names(bare, "combined")
