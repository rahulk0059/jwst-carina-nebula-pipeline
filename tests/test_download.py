"""Tests for the MAST downloader, all staying off the network."""
from __future__ import annotations

import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np
import pytest
from astropy.table import Table

import csv

from jwst_stack import download


def _products_table(rows):
    return Table(
        {
            "productFilename": [r[0] for r in rows],
            "productSubGroupDescription": [r[1] for r in rows],
            "filters": [r[2] for r in rows],
            "dataURI": [r[3] for r in rows],
            "size": [r[4] for r in rows],
            "obs_id": [r[5] for r in rows],
        }
    )


CAL_ROWS = [
    ("jw02731001001_02105_00001_nrca1_cal.fits", "CAL", "F200W", "mast:JWST/product/a1", 117573120, "jw02731-o001_t017_nircam_clear-f200w"),
    ("jw02731001003_02105_00002_nrcblong_cal.fits", "CAL", "F444W;F470N", "mast:JWST/product/b1", 62000, "jw02731-o001_t017_nircam_f444w-f470n"),
    ("jw02731001002_02105_00004_nrcb4_i2d.fits", "I2D", "F200W", "mast:JWST/product/c1", 118702080, "jw02731-o001_t017_nircam_clear-f200w"),
    ("jw02731001001_02105_00005_nrcalong_i2d.fits", "I2D", "F444W;F470N", "mast:JWST/product/d1", 100, "jw02731-o001_t017_nircam_f444w-f470n"),
]


def test_parse_product_name_sw_and_lw():
    assert download.parse_product_name("jw02731001001_02105_00001_nrca1_cal.fits") == ("nrca1", "cal")
    assert download.parse_product_name("jw02731001001_02105_00001_nrca4_i2d.fits") == ("nrca4", "i2d")
    assert download.parse_product_name("jw02731001003_02105_00002_nrcblong_cal.fits") == ("nrcblong", "cal")
    assert download.parse_product_name("jw02731001002_02105_00004_nrcalong_i2d.fits") == ("nrcalong", "i2d")
    assert download.parse_product_name("jw02731001001_02105_00001_preview.jpg") is None


def test_split_filter_pupil():
    assert download.split_filter_pupil("F200W") == ("F200W", "CLEAR")
    assert download.split_filter_pupil("F444W;F470N") == ("F444W", "F470N")
    assert download.split_filter_pupil(" F090W; CLEAR ") == ("F090W", "CLEAR")


def test_visit_and_exposure_from_name():
    assert download.visit_and_exposure_from_name("jw02731001003_02105_00002_nrca1_cal.fits") == (3, 2)
    assert download.visit_and_exposure_from_name("jw02731001001_02105_00005_nrcalong_cal.fits") == (1, 5)


def test_build_items_selection_and_sort():
    prods = _products_table(CAL_ROWS)
    items = download.build_items(prods, {"CAL"}, "/tmp/dst")
    assert [i.filename for i in items] == [
        "jw02731001001_02105_00001_nrca1_cal.fits",
        "jw02731001003_02105_00002_nrcblong_cal.fits",
    ]
    assert items[0].filter_name == "F200W"
    assert items[0].pupil == "CLEAR"
    assert items[0].visit == 1
    assert items[0].exposure == 1
    assert items[0].detector == "nrca1"
    assert items[0].expected_size == 117573120
    assert items[1].filter_name == "F444W"
    assert items[1].pupil == "F470N"
    assert items[1].url.startswith("https://mast.stsci.edu/api/v0.1/Download/file?uri=")


def test_build_items_filters_detectors_kinds():
    prods = _products_table(CAL_ROWS)
    items = download.build_items(prods, {"CAL", "I2D"}, "/tmp/dst", filter_names={"F200W"})
    assert {i.kind for i in items} == {"CAL", "I2D"}
    assert all(i.filter_name == "F200W" for i in items)
    items = download.build_items(prods, {"CAL"}, "/tmp/dst", detectors={"nrca1"})
    assert [i.filename for i in items] == ["jw02731001001_02105_00001_nrca1_cal.fits"]


def test_check_disk_too_small(monkeypatch, tmp_path):
    class FakeUsage:
        free = 1000

    monkeypatch.setattr(download.shutil, "disk_usage", lambda p: FakeUsage())
    free, ok, msg = download.check_disk(tmp_path, total_bytes=5000, safety=1.1)
    assert not ok
    assert "GB free" in msg
    assert free == 1000


def test_check_disk_enough(monkeypatch, tmp_path):
    class FakeUsage:
        free = 10_000

    monkeypatch.setattr(download.shutil, "disk_usage", lambda p: FakeUsage())
    _, ok, _ = download.check_disk(tmp_path, total_bytes=100, safety=2.0)
    assert ok


def test_file_size_matches(tmp_path):
    p = tmp_path / "f"
    assert not download.file_size_matches(p, 5)
    p.write_bytes(b"12345")
    assert download.file_size_matches(p, 5)
    assert not download.file_size_matches(p, 4)
    assert not download.file_size_matches(p, 6)


class _RangeHandler(BaseHTTPRequestHandler):
    payload = b""

    def do_GET(self):
        start = 0
        rng = self.headers.get("Range")
        if rng:
            start = int(rng.split("=")[1].split("-")[0])
        data = self.payload[start:]
        if start > 0:
            self.send_response(206, "Partial Content")
        else:
            self.send_response(200, "OK")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def log_message(self, *args):
        pass


@pytest.fixture
def range_server():
    payload = (bytes(range(256)) * 2000)[:256000]
    _RangeHandler.payload = payload
    server = ThreadingHTTPServer(("127.0.0.1", 0), _RangeHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    url = f"http://127.0.0.1:{server.server_address[1]}/file"
    yield url, payload
    server.shutdown()


def test_download_url_skip_existing(tmp_path, range_server):
    url, payload = range_server
    dest = tmp_path / "f"
    dest.write_bytes(payload)
    status = download.download_url(url, dest, len(payload), resume=True)
    assert status == "skipped"


def test_download_url_resume_from_partial(tmp_path, range_server):
    url, payload = range_server
    dest = tmp_path / "f"
    dest.write_bytes(payload[:1000])
    status = download.download_url(url, dest, len(payload), resume=True)
    assert status == "resumed"
    assert dest.read_bytes() == payload


def test_download_url_restarts_when_existing_larger(tmp_path, range_server):
    url, payload = range_server
    dest = tmp_path / "f"
    dest.write_bytes(payload + b"junk")
    status = download.download_url(url, dest, len(payload), resume=True)
    assert status in ("ok", "resumed")
    assert dest.read_bytes() == payload


class _IgnoresRangeHandler(BaseHTTPRequestHandler):
    """Answers every request with 200 and the full body, ignoring Range."""

    payload = b""

    def do_GET(self):
        data = self.payload
        self.send_response(200, "OK")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def log_message(self, *args):
        pass


@pytest.fixture
def range_ignoring_server():
    payload = (bytes(range(256)) * 2000)[:256000]
    _IgnoresRangeHandler.payload = payload
    server = ThreadingHTTPServer(("127.0.0.1", 0), _IgnoresRangeHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    url = f"http://127.0.0.1:{server.server_address[1]}/file"
    yield url, payload
    server.shutdown()


def test_download_url_truncates_when_server_ignores_range(tmp_path, range_ignoring_server):
    """A 200 answer to a Range request must restart, not append.

    Opening the partial in "ab" and calling seek(0) is a no-op for writes
    (O_APPEND), so the full body used to be appended onto the partial and left
    an oversized corrupt file behind.
    """
    url, payload = range_ignoring_server
    dest = tmp_path / "f"
    dest.write_bytes(payload[:1000])

    status = download.download_url(url, dest, len(payload), resume=True)

    assert status == "ok"
    assert dest.stat().st_size == len(payload)
    assert dest.read_bytes() == payload


def test_download_url_truncates_when_server_ignores_range_no_resume(
    tmp_path, range_ignoring_server
):
    url, payload = range_ignoring_server
    dest = tmp_path / "f"
    dest.write_bytes(payload[:1000])

    status = download.download_url(url, dest, len(payload), resume=False)

    assert status == "ok"
    assert dest.read_bytes() == payload


def test_download_items_writes_manifest(tmp_path, monkeypatch):
    item = download.DownloadItem(
        filename="x.fits",
        url="http://localhost/x",
        expected_size=10,
        filter_name="F200W",
        pupil="CLEAR",
        detector="nrca1",
        visit=1,
        exposure=1,
        kind="CAL",
        obs_id="jw-ob",
    )

    def fake_download(file_url, dest, expected, resume=True, **kw):
        dest.write_bytes(b"0" * expected)
        return "ok"

    monkeypatch.setattr(download, "download_url", fake_download)
    dest = tmp_path / "data"
    manifest = tmp_path / "manifest.csv"
    items = download.download_items([item], dest, manifest)
    assert items[0].status == "ok"
    assert (dest / "x.fits").exists()
    rows = download.read_manifest(manifest)
    assert len(rows) == 1
    assert rows[0]["filter_name"] == "F200W"
    assert rows[0]["pupil"] == "CLEAR"
    assert rows[0]["detector"] == "nrca1"
    assert rows[0]["visit"] == "1"
    assert rows[0]["status"] == "ok"
    assert rows[0]["local_path"] == str(dest / "x.fits")


def test_run_download_plan_only_never_prompts_or_downloads(tmp_path, monkeypatch, capsys):
    """--plan-only must print the plan and stop, not fall through to the prompt."""
    obs = Table({"obs_id": ["jw02731-o001_t017_nircam_clear-f200w"]})
    prods = _products_table(CAL_ROWS[:1])

    monkeypatch.setattr(download, "query_obs_table", lambda **kw: obs)
    monkeypatch.setattr(download, "query_products", lambda rows: prods)
    monkeypatch.setattr(download.shutil, "disk_usage", lambda p: type("U", (), {"free": 1 << 40})())

    def explode(*a, **k):  # pragma: no cover - must not be reached
        raise AssertionError("plan-only attempted a download")

    monkeypatch.setattr(download, "download_url", explode)

    import argparse

    args = argparse.Namespace(
        input_dir=str(tmp_path / "data"),
        outdir=str(tmp_path / "out"),
        manifest=str(tmp_path / "manifest.csv"),
        filters=[],
        detectors=[],
        products="cal,i2d",
        proposal_id=2731,
        yes=False,
        plan_only=True,
    )

    def never(msg):  # pragma: no cover - must not be reached
        raise AssertionError("plan-only prompted for confirmation")

    assert download.run_download(args, confirm=never) == 0
    out = capsys.readouterr().out
    assert "Planned downloads:" in out
    assert "plan-only: nothing downloaded" in out
    data = tmp_path / "data"
    assert not data.exists() or not any(data.rglob("*"))


def test_run_download_cancel_and_proceed(tmp_path, monkeypatch):
    obs = Table({"obs_id": ["jw02731-o001_t017_nircam_clear-f200w"]})
    prods = _products_table(CAL_ROWS[:1])

    monkeypatch.setattr(download, "query_obs_table", lambda **kw: obs)
    monkeypatch.setattr(download, "query_products", lambda rows: prods)
    monkeypatch.setattr(download.shutil, "disk_usage", lambda p: type("U", (), {"free": 1 << 40})())

    def fake_download(file_url, dest, expected, resume=True, **kw):
        dest.write_bytes(b"0" * expected)
        return "ok"

    monkeypatch.setattr(download, "download_url", fake_download)

    import argparse

    def make_args(yes):
        return argparse.Namespace(
            input_dir=str(tmp_path / "data"),
            outdir=str(tmp_path / "out"),
            manifest=str(tmp_path / "manifest.csv"),
            filters=[],
            detectors=[],
            products="cal,i2d",
            proposal_id=2731,
            yes=yes,
            plan_only=False,
        )

    consensus = []
    base = make_args(False)
    code = download.run_download(base, confirm=lambda msg: consensus.append(msg) or "n")
    assert code == 0
    assert not (tmp_path / "data" / "jw02731001001_02105_00001_nrca1_cal.fits").exists()
    assert consensus, "confirmation was prompted"

    code = download.run_download(make_args(True))
    assert code == 0
    assert (tmp_path / "data" / "jw02731001001_02105_00001_nrca1" / "jw02731001001_02105_00001_nrca1_cal.fits").exists()
    assert (tmp_path / "manifest.csv").exists()


def test_summarize_items_includes_total():
    prods = _products_table(CAL_ROWS)
    items = download.build_items(prods, {"CAL", "I2D"}, "/tmp/dst")
    out = download.summarize_items(items)
    assert "Total" in out
    assert "F200W" in out


def _make_row(filename, subgroup, filters="F200W", size=1000):
    return (filename, subgroup, filters, f"mast:JWST/product/{filename}", size, "jw-ob")


def test_is_combined_i2d():
    assert download.is_combined_i2d("jw02731-o001_t017_nircam_clear-f200w_i2d.fits")
    assert not download.is_combined_i2d("jw02731001001_02105_00001_nrca1_i2d.fits")


def test_plan_stage2_selection_and_existing(tmp_path):
    rows = [
        _make_row("jw02731001008_02105_00004_nrca1_cal.fits", "CAL", size=100),  # visit 8 cal
        _make_row("jw02731001002_02105_00003_nrca3_cal.fits", "CAL", size=100),
        _make_row("jw02731-o001_t017_nircam_clear-f200w_i2d.fits", "I2D", size=200),
        _make_row("jw02731001001_02105_00005_nrca1_i2d.fits", "I2D", size=110),  # visit 1 subset
        _make_row("jw02731001001_02105_00006_nrcb4_i2d.fits", "I2D", size=110),  # visit 1 subset
        _make_row("jw02731001003_02105_00001_nrca1_i2d.fits", "I2D", size=110),  # visit 3 nrca1 subset
        _make_row("jw02731001002_02105_00005_nrca2_i2d.fits", "I2D", size=110),  # visit 2 nrca2 NOT subset
    ]
    prods = _products_table(rows)
    data_root = tmp_path / "data"
    i2d_root = tmp_path / "i2d"
    data_root.mkdir()
    i2d_root.mkdir()
    # Pre-seed one subset i2d at the exact size -> counts as existing.
    (i2d_root / "jw02731001001_02105_00005_nrca1_i2d.fits").write_bytes(b"x" * 110)
    # Pre-seed one cal at the exact size -> counts as existing.
    (data_root / "jw02731001008_02105_00004_nrca1_cal.fits").write_bytes(b"x" * 100)
    # Pre-seed one cal at the wrong size -> must stay 'new'.
    (data_root / "jw02731001002_02105_00003_nrca3_cal.fits").write_bytes(b"x" * 1)

    plan = download.plan_stage2(prods, data_root, i2d_root)
    assert len(plan["all_cal"]) == 2
    assert len(plan["all_i2d"]) == 5
    assert [i.filename for i in plan["combined"]] == ["jw02731-o001_t017_nircam_clear-f200w_i2d.fits"]
    assert sorted(i.filename for i in plan["subset"]) == [
        "jw02731001001_02105_00005_nrca1_i2d.fits",
        "jw02731001001_02105_00006_nrcb4_i2d.fits",
        "jw02731001003_02105_00001_nrca1_i2d.fits",
    ]

    cal_existing = {it.filename for it, s in plan["categories"]["cal"] if s == "existing"}
    cal_new = {it.filename for it, s in plan["categories"]["cal"] if s == "new"}
    assert "jw02731001002_02105_00003_nrca3_cal.fits" not in cal_existing
    assert "jw02731001002_02105_00003_nrca3_cal.fits" in cal_new
    assert "jw02731001008_02105_00004_nrca1_cal.fits" in cal_existing

    sub_existing = [it.filename for it, s in plan["categories"]["subset_i2d"] if s == "existing"]
    assert sub_existing == ["jw02731001001_02105_00005_nrca1_i2d.fits"]


def test_format_stage2_plan_reports_scratch():
    plan = {"categories": {"cal": [], "combined_i2d": [], "subset_i2d": []}}
    text = download.format_stage2_plan(plan, scratch_estimate_gb=9.5)
    assert "est. mosaic scratch: 9.5 GB" in text
    assert "to download" in text


def _two_filter_stage2_table():
    """Products for two filters, one combined i2d each, plus a shared subset."""
    return _products_table(
        [
            _make_row("jw02731001001_02105_00001_nrca1_cal.fits", "CAL", "F200W", size=100),
            _make_row("jw02731001001_02105_00001_nrca1_cal.fits", "CAL", "F187N", size=100),
            _make_row("jw02731-o001_t017_nircam_clear-f200w_i2d.fits", "I2D", "F200W", size=200),
            _make_row("jw02731-o001_t017_nircam_clear-f187n_i2d.fits", "I2D", "F187N", size=200),
            _make_row("jw02731001001_02105_00005_nrca1_i2d.fits", "I2D", "F187N", size=110),
        ]
    )


def test_plan_stage2_respects_filters(tmp_path):
    """--stage2 must plan the requested filter, not a hard-coded F200W."""
    prods = _two_filter_stage2_table()
    data_root = tmp_path / "data"
    i2d_root = tmp_path / "i2d"
    data_root.mkdir()
    i2d_root.mkdir()

    f187n = download.plan_stage2(prods, data_root, i2d_root, ["F187N"])
    assert f187n["filters"] == ["F187N"]
    assert [i.filter_name for i in f187n["all_cal"]] == ["F187N"]
    assert [i.filename for i in f187n["combined"]] == [
        "jw02731-o001_t017_nircam_clear-f187n_i2d.fits"
    ]
    assert "F200W" not in download.format_stage2_plan(f187n)

    # the default is only a fallback for an unspecified filter
    default = download.plan_stage2(prods, data_root, i2d_root)
    assert default["filters"] == ["F200W"]
    assert [i.filter_name for i in default["all_cal"]] == ["F200W"]

    # both filters at once merge, and the label says so
    both = download.plan_stage2(prods, data_root, i2d_root, ["F187N", "F200W"])
    assert both["filters"] == ["F187N", "F200W"]
    assert sorted(i.filter_name for i in both["all_cal"]) == ["F187N", "F200W"]
    text = download.format_stage2_plan(both)
    assert "F187N,F200W cal files (all)" in text
    assert "F200W cal files (all 160)" not in text


def test_parse_stage2_kinds_validates_and_orders():
    assert download._parse_stage2_kinds("cal,combined_i2d") == ["cal", "combined_i2d"]
    # order follows the canonical order, not the order given
    assert download._parse_stage2_kinds("subset_i2d,cal") == ["cal", "subset_i2d"]
    assert download._parse_stage2_kinds(" cal , subset_i2d ") == ["cal", "subset_i2d"]
    with pytest.raises(ValueError, match="unknown --stage2-kinds"):
        download._parse_stage2_kinds("cal,bogus")
    with pytest.raises(ValueError, match="at least one"):
        download._parse_stage2_kinds(",")


def test_stage2_kinds_excludes_the_perexposure_subset(monkeypatch, capsys, tmp_path):
    """cal + combined_i2d only: the 161-file F187N scope, via the CLI."""
    monkeypatch.setattr(download, "query_obs_table", lambda **kw: _obs_table())
    monkeypatch.setattr(download, "query_products", lambda rows: _two_filter_stage2_table())
    monkeypatch.setattr(download, "check_disk", lambda *a, **k: (1e12, True, "ok"))

    args = download.parse_args(
        [
            "--stage2",
            "--plan-only",
            "--filters",
            "F187N",
            "--stage2-kinds",
            "cal,combined_i2d",
            "--input-dir",
            str(tmp_path / "data"),
            "--i2d-dir",
            str(tmp_path / "i2d"),
            "--outdir",
            str(tmp_path),
        ]
    )
    assert download.run_stage2(args) == 0
    out = capsys.readouterr().out
    assert "F187N cal files (all)" in out
    assert "F187N combined i2d mosaic" in out
    assert "validation subset" not in out
    assert "F200W" not in out


def test_download_url_reports_incremental_progress(tmp_path, range_server):
    """A large transfer must be observable, not silent until it finishes."""
    url, payload = range_server
    dest = tmp_path / "big.fits"
    seen = []
    status = download.download_url(
        url, dest, len(payload), resume=True,
        progress=lambda w, t: seen.append((w, t)), report_every=16,
    )
    assert status == "ok"
    assert seen, "progress callback was never invoked"
    assert all(w <= len(payload) for w, _ in seen)
    assert seen == sorted(seen), "progress went backwards"
    assert seen[-1][0] == len(payload)


def test_download_url_abandons_a_stalled_transfer(tmp_path, range_server):
    """A transfer that never finishes must be reported, not run forever."""
    url, payload = range_server
    dest = tmp_path / "stalled.fits"
    status = download.download_url(
        url, dest, len(payload) * 4, resume=True, max_seconds=-1.0
    )
    assert status.startswith("error:abandoned")
    assert dest.stat().st_size < len(payload) * 4


def test_item_failed_treats_only_ok_and_skipped_as_success():
    def mk(status):
        return download.DownloadItem(
            filename="f.fits", kind="CAL", filter_name="F200W", pupil="CLEAR",
            visit=1, detector="nrca1", exposure=1, obs_id="o1", url="u",
            expected_size=10, status=status,
        )

    assert not download.item_failed(mk("ok"))
    assert not download.item_failed(mk("skipped"))
    for bad in ("error:<urlopen error [Errno 11001] getaddrinfo failed>",
                "error:The read operation timed out", "size_mismatch", ""):
        assert download.item_failed(mk(bad)), bad


def test_run_stage2_returns_nonzero_when_a_file_fails(monkeypatch, capsys, tmp_path):
    """A partial download must not report itself complete with exit 0."""
    monkeypatch.setattr(download, "query_obs_table", lambda **kw: _obs_table())
    monkeypatch.setattr(download, "query_products", lambda rows: _two_filter_stage2_table())
    monkeypatch.setattr(download, "check_disk", lambda *a, **k: (1e12, True, "ok"))
    monkeypatch.setattr(download, "_load_grid_shape", lambda outdir: (100, 100))

    calls = {"n": 0}
    real = download.download_items

    def flaky(items, dest_root, manifest_path, progress=print):
        calls["n"] += 1
        out = real(items, dest_root, manifest_path, progress)
        if calls["n"] == 1:
            out[0].status = "error:<urlopen error [Errno 11001] getaddrinfo failed>"
        return out

    monkeypatch.setattr(download, "download_url",
                        lambda url, dest, size, resume=False, **kw: "ok")
    monkeypatch.setattr(download, "download_items", flaky)

    args = download.parse_args(
        [
            "--stage2", "--yes", "--filters", "F187N",
            "--stage2-kinds", "cal,combined_i2d",
            "--input-dir", str(tmp_path / "data"),
            "--i2d-dir", str(tmp_path / "i2d"),
            "--outdir", str(tmp_path),
        ]
    )
    assert download.run_stage2(args) == 1
    out = capsys.readouterr().out
    assert "INCOMPLETE" in out
    assert "Re-run the same command to resume" in out
    assert "Stage-2 download complete" not in out


def test_run_download_returns_nonzero_when_a_file_fails(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(download, "query_obs_table", lambda **kw: _obs_table())
    monkeypatch.setattr(download, "query_products", lambda rows: _products_table(CAL_ROWS))
    monkeypatch.setattr(download, "check_disk", lambda *a, **k: (1e12, True, "ok"))
    monkeypatch.setattr(download, "download_url",
                        lambda url, dest, size, resume=False, **kw: "error:boom")

    args = download.parse_args(
        ["--yes", "--filters", "F200W", "--products", "cal",
         "--input-dir", str(tmp_path / "data"), "--outdir", str(tmp_path)]
    )
    assert download.run_download(args) == 1
    assert "INCOMPLETE" in capsys.readouterr().out


def test_resolve_stage2_filters_prefers_requested():
    def make(filters):
        return download.argparse.Namespace(filters=filters)

    assert download.resolve_stage2_filters(make(["f187n"])) == ["F187N"]
    assert download.resolve_stage2_filters(make(["F187N", "F090W"])) == ["F187N", "F090W"]
    assert download.resolve_stage2_filters(make([])) == ["F200W"]
    assert download.resolve_stage2_filters(make(None)) == ["F200W"]


def test_run_stage2_plan_output_labels_requested_filter(monkeypatch, capsys, tmp_path):
    """End-to-end: the printed plan says F187N, never F200W."""
    monkeypatch.setattr(download, "query_obs_table", lambda **kw: _obs_table())
    monkeypatch.setattr(download, "query_products", lambda rows: _two_filter_stage2_table())
    monkeypatch.setattr(download, "check_disk", lambda *a, **k: (1e12, True, "ok"))

    args = download.parse_args(
        [
            "--stage2",
            "--plan-only",
            "--filters",
            "F187N",
            "--input-dir",
            str(tmp_path / "data"),
            "--i2d-dir",
            str(tmp_path / "i2d"),
            "--outdir",
            str(tmp_path),
        ]
    )
    assert download.run_stage2(args) == 0
    out = capsys.readouterr().out
    assert "F187N cal files (all)" in out
    assert "F187N combined i2d mosaic" in out
    assert "jw02731-o001_t017_nircam_clear-f187n_i2d.fits" not in out
    assert "F200W" not in out


def test_estimate_mosaic_scratch_gb():
    gb = download.estimate_mosaic_scratch_gb((100, 100), 4, overlay_naxis=100)
    assert gb == pytest.approx((100 * 100 * 10 + 4 * 100 * 100 * 10) / 1e9)


def test_stage2_main_dispatch(monkeypatch):
    called = {}

    def fake_stage2(a, confirm=None):
        called["ran"] = True
        return 42

    monkeypatch.setattr(download, "run_stage2", fake_stage2)
    code = download.main(["--stage2", "--plan-only"])
    assert code == 42
    assert called["ran"]


VERIFY_ROWS = [
    ("jw02731001001_02105_00001_nrca1_cal.fits", "CAL", "F200W", "mast:JWST/product/v1", 1000, "obs-a"),
    ("jw02731001001_02105_00002_nrca1_cal.fits", "CAL", "F200W", "mast:JWST/product/v2", 2000, "obs-a"),
    ("jw02731001001_02105_00003_nrca1_cal.fits", "CAL", "F200W", "mast:JWST/product/v3", 3000, "obs-a"),
    ("jw02731001001_02105_00001_nrca1_i2d.fits", "I2D", "F200W", "mast:JWST/product/v4", 400, "obs-a"),
]


def _obs_table():
    return Table({"obs_id": ["jw02731-o001_t017_nircam_clear-f200w"]})


def _place(root, filename, nbytes):
    group = download._group_dir_name(filename)
    path = root / group / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\0" * nbytes)
    return path


MULTI_FILTER_ROWS = [
    ("jw02731001001_02105_00001_nrca1_cal.fits", "CAL", "F200W", "mast:JWST/product/a1", 117573120, "jw02731-o001_t017_nircam_clear-f200w"),
    ("jw02731001001_02103_00001_nrca1_cal.fits", "CAL", "F090W", "mast:JWST/product/a2", 117573120, "jw02731-o001_t017_nircam_clear-f090w"),
    ("jw02731-o001_t017_nircam_clear-f090w_i2d.fits", "I2D", "F090W", "mast:JWST/product/c2", 5415445, "jw02731-o001_t017_nircam_clear-f090w"),
]


def test_verify_covers_every_filter_present_on_disk(tmp_path):
    """A plan missing a filter must not report that filter's files as 'extra'.

    The data root holds more than one filter, so a single-filter plan labels
    every other filter's files 'extra' and the non-zero exit stops gating
    anything real.
    """
    prods = _products_table(MULTI_FILTER_ROWS)
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    for name, size in (
        ("jw02731001001_02105_00001_nrca1_cal.fits", 117573120),
        ("jw02731001001_02103_00001_nrca1_cal.fits", 117573120),
    ):
        p = cal_root / "grp" / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"0" * size)
    ip = i2d_root / "grp" / "jw02731-o001_t017_nircam_clear-f090w_i2d.fits"
    ip.parent.mkdir(parents=True, exist_ok=True)
    ip.write_bytes(b"0" * (5415445 // 2))  # size mismatch on purpose

    # Both filters: F090W is recognised, and its wrong size is caught.
    cal, i2d = download.verify_products(
        prods, cal_root, i2d_root, filter_name=["F200W", "F090W"]
    )
    statuses = {r.item.filename: r.status for r in cal.rows}
    assert statuses["jw02731001001_02105_00001_nrca1_cal.fits"] == "ok"
    assert statuses["jw02731001001_02103_00001_nrca1_cal.fits"] == "ok"
    assert "extra" not in statuses.values()
    istatus = {r.item.filename: r.status for r in i2d.rows}
    assert istatus["jw02731-o001_t017_nircam_clear-f090w_i2d.fits"] == "size_mismatch"

    # A single-filter plan is the trap: the other filter becomes 'extra'.
    cal1, _ = download.verify_products(
        prods, cal_root, i2d_root, filter_name=["F200W"]
    )
    statuses1 = {r.item.filename: r.status for r in cal1.rows}
    assert statuses1["jw02731001001_02103_00001_nrca1_cal.fits"] == "extra"


def test_nircam_filters_lists_every_filter_in_the_plan():
    prods = _products_table(MULTI_FILTER_ROWS)
    assert download.nircam_filters(prods) == ["F090W", "F200W"]


def test_verify_includes_combined_mosaic(tmp_path):
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02105_00001_nrca1_cal.fits", 1000)  # ok
    _place(cal_root, "jw02731001001_02105_00002_nrca1_cal.fits", 1500)  # truncated
    # ...00003 deliberately absent
    _place(cal_root, "jw02731001999_02105_00009_nrca1_cal.fits", 77)  # extra
    _place(i2d_root, "jw02731001001_02105_00001_nrca1_i2d.fits", 400)  # ok

    cal, i2d = download.verify_products(
        _products_table(VERIFY_ROWS), cal_root, i2d_root
    )

    assert (cal.count("ok"), cal.count("size_mismatch")) == (1, 1)
    assert (cal.count("missing"), cal.count("extra")) == (1, 1)
    assert i2d.count("ok") == 1 and not i2d.problems

    bad = next(r for r in cal.rows if r.status == "size_mismatch")
    assert bad.item.expected_size == 2000 and bad.actual_size == 1500
    assert bad.item.filename == "jw02731001001_02105_00002_nrca1_cal.fits"


def test_verify_report_text_lists_problems(tmp_path):
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02105_00002_nrca1_cal.fits", 10)
    cal, i2d = download.verify_products(
        _products_table(VERIFY_ROWS), cal_root, i2d_root
    )
    text = download.format_verify_report([cal, i2d])
    assert "size_mismatch" in text
    assert "missing" in text
    assert "VERIFY FAILED" in text
    assert "1,500" not in text  # the truncated one is 10 bytes, not 1500
    assert "jw02731001001_02105_00002_nrca1_cal.fits" in text


def test_write_verify_manifest_records_actual_sizes(tmp_path):
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02105_00002_nrca1_cal.fits", 1500)
    cal, _ = download.verify_products(_products_table(VERIFY_ROWS), cal_root, i2d_root)

    path = download.write_verify_manifest(cal, tmp_path / "m.csv")
    rows = {r["filename"]: r for r in download.read_manifest(path)}

    assert rows["jw02731001001_02105_00001_nrca1_cal.fits"]["status"] == "missing"
    assert rows["jw02731001001_02105_00001_nrca1_cal.fits"]["size"] == "0"
    assert rows["jw02731001001_02105_00002_nrca1_cal.fits"]["size"] == "1500"
    assert rows["jw02731001001_02105_00002_nrca1_cal.fits"]["status"] == "size_mismatch"


def test_verify_includes_combined_mosaic(tmp_path):
    """The combined i2d mosaic has no detector/exposure in its name.

    ``build_items`` drops such products, which made verify report the mosaic
    as an unexplained "extra" instead of comparing it with the MAST size.
    """
    combined = "jw02731-o001_t017_nircam_clear-f200w_i2d.fits"
    rows = VERIFY_ROWS + [
        (combined, "I2D", "F200W", "mast:JWST/product/v5", 5000, "obs-a")
    ]
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(i2d_root, combined, 2000)

    cal, i2d = download.verify_products(
        _products_table(rows), cal_root, i2d_root
    )

    assert i2d.count("extra") == 0
    row = next(r for r in i2d.rows if r.item.filename == combined)
    assert row.status == "size_mismatch"
    assert row.item.expected_size == 5000 and row.actual_size == 2000
    assert row.item.detector == "combined"


def test_run_verify_rebuilds_manifests_and_returns_failure(tmp_path, monkeypatch):
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02105_00001_nrca1_cal.fits", 1000)
    outdir = tmp_path / "out"
    table = _products_table(VERIFY_ROWS)

    monkeypatch.setattr(
        download, "query_obs_table", lambda proposal_id: _obs_table()
    )
    monkeypatch.setattr(download, "query_products", lambda rows: table)
    monkeypatch.setattr(
        download, "is_nircam_obs", lambda o: True
    )

    args = download.parse_args(
        [
            "--verify",
            "--input-dir",
            str(cal_root),
            "--i2d-dir",
            str(i2d_root),
            "--outdir",
            str(outdir),
        ]
    )
    code = download.run_verify(args)

    assert code == 1  # there are missing files
    assert (outdir / "verify_report.txt").exists()
    cal_rows = download.read_manifest(outdir / "stage2_cal_manifest.csv")
    assert {r["filename"] for r in cal_rows} == {r[0] for r in VERIFY_ROWS if r[1] == "CAL"}
    assert next(r for r in cal_rows if r["status"] == "ok")["size"] == "1000"
    i2d_rows = download.read_manifest(outdir / "stage2_i2d_manifest.csv")
    assert {r["filename"] for r in i2d_rows} == {r[0] for r in VERIFY_ROWS if r[1] == "I2D"}


def test_run_verify_no_rebuild_leaves_manifests(tmp_path, monkeypatch):
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    outdir = tmp_path / "out"
    monkeypatch.setattr(
        download, "query_obs_table", lambda proposal_id: _obs_table()
    )
    monkeypatch.setattr(download, "query_products", lambda rows: _products_table(VERIFY_ROWS))
    monkeypatch.setattr(download, "is_nircam_obs", lambda o: True)

    args = download.parse_args(
        [
            "--verify",
            "--no-rebuild",
            "--input-dir",
            str(cal_root),
            "--i2d-dir",
            str(i2d_root),
            "--outdir",
            str(outdir),
        ]
    )
    assert download.run_verify(args) == 1
    assert (outdir / "verify_report.txt").exists()
    assert not (outdir / "stage2_cal_manifest.csv").exists()


def test_verify_main_dispatch(monkeypatch):
    monkeypatch.setattr(download, "run_verify", lambda a: 7)
    assert download.main(["--verify"]) == 7


# --------------------------------------------------------------------------
# --require: scope the verify exit code to one filter without narrowing the scan
# --------------------------------------------------------------------------

SCOPE_ROWS = [
    # F090W cal, present with the exact size -> ok
    ("jw02731001001_02103_00001_nrca1_cal.fits", "CAL", "F090W", "mast:JWST/product/a1", 1000, "obs-f090w"),
    # F090W cal, absent -> must gate when F090W is the required filter
    ("jw02731001001_02103_00002_nrca1_cal.fits", "CAL", "F090W", "mast:JWST/product/a2", 2000, "obs-f090w"),
    # F200W cal, absent -> must NOT gate when only F090W is required
    ("jw02731001001_02105_00001_nrca1_cal.fits", "CAL", "F200W", "mast:JWST/product/b1", 3000, "obs-f200w"),
    # F200W cal, present but truncated -> a corrupt file gates at any scope
    ("jw02731001001_02105_00002_nrca1_cal.fits", "CAL", "F200W", "mast:JWST/product/b2", 4000, "obs-f200w"),
]


def _scope_fixture(tmp_path):
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02103_00001_nrca1_cal.fits", 1000)  # ok
    _place(cal_root, "jw02731001001_02105_00002_nrca1_cal.fits", 1500)  # truncated
    return download.verify_products(
        _products_table(SCOPE_ROWS), cal_root, i2d_root,
        filter_name=["F090W", "F200W"],
    )


def test_verify_scope_gates_absent_products_only_for_named_filters(tmp_path):
    cal, _ = _scope_fixture(tmp_path)
    gating, relaxed = cal.split_by_scope({"F090W"})
    gated = {(r.status, r.item.filter_name) for r in gating}
    assert ("missing", "F090W") in gated
    assert ("size_mismatch", "F200W") in gated
    assert ("missing", "F200W") in {(r.status, r.item.filter_name) for r in relaxed}


def test_verify_scope_never_relaxes_a_corrupt_or_unexpected_file(tmp_path):
    """size_mismatch / extra gate at any scope: they are disk anomalies."""
    cal, _ = _scope_fixture(tmp_path)
    for scope in ({"F090W"}, {"F200W"}, None):
        gating, _relaxed = cal.split_by_scope(scope)
        assert any(r.status == "size_mismatch" for r in gating), scope


def test_verify_without_scope_gates_every_problem(tmp_path):
    """The default whole-program behaviour is unchanged."""
    cal, _ = _scope_fixture(tmp_path)
    gating, relaxed = cal.split_by_scope(None)
    assert len(gating) == len(cal.problems)
    assert relaxed == []


def test_verify_report_text_separates_out_of_scope_rows(tmp_path):
    cal, i2d = _scope_fixture(tmp_path)
    text = download.format_verify_report([cal, i2d], scope={"F090W"})
    assert "kinds" not in text.split("problems (gating)")[0]
    assert "gating scope: filters F090W" in text
    assert "outside the gating scope" in text
    # the relaxed F200W absence is still shown, just not counted
    assert "jw02731001001_02105_00001_nrca1_cal.fits" in text
    assert "do not affect the exit code" in text


def _scope_args(cal_root, i2d_root, outdir, require=None):
    argv = [
        "--verify",
        "--input-dir", str(cal_root),
        "--i2d-dir", str(i2d_root),
        "--outdir", str(outdir),
    ]
    if require:
        argv += ["--require", *require]
    return download.parse_args(argv)


def _patch_mast(monkeypatch, rows):
    monkeypatch.setattr(download, "query_obs_table", lambda proposal_id: _obs_table())
    monkeypatch.setattr(download, "query_products", lambda r: _products_table(rows))
    monkeypatch.setattr(download, "is_nircam_obs", lambda o: True)


def test_run_verify_require_exits_zero_for_a_clean_filter(tmp_path, monkeypatch):
    """A fully-downloaded filter exits 0 even though other filters are absent."""
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02103_00001_nrca1_cal.fits", 1000)
    outdir = tmp_path / "out"
    _patch_mast(
        monkeypatch,
        [
            ("jw02731001001_02103_00001_nrca1_cal.fits", "CAL", "F090W", "mast:JWST/product/a1", 1000, "obs-f090w"),
            ("jw02731001001_02105_00001_nrca1_cal.fits", "CAL", "F200W", "mast:JWST/product/b1", 3000, "obs-f200w"),
        ],
    )

    assert download.run_verify(_scope_args(cal_root, i2d_root, outdir, ["F090W"])) == 0
    # ... and the unscoped run still fails, because F200W is absent.
    assert download.run_verify(_scope_args(cal_root, i2d_root, outdir)) == 1


def test_run_verify_require_still_fails_on_an_incomplete_filter(tmp_path, monkeypatch):
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02103_00001_nrca1_cal.fits", 1000)  # 2nd F090W absent
    outdir = tmp_path / "out"
    _patch_mast(
        monkeypatch,
        [
            ("jw02731001001_02103_00001_nrca1_cal.fits", "CAL", "F090W", "mast:JWST/product/a1", 1000, "obs-f090w"),
            ("jw02731001001_02103_00002_nrca1_cal.fits", "CAL", "F090W", "mast:JWST/product/a2", 2000, "obs-f090w"),
        ],
    )

    assert download.run_verify(_scope_args(cal_root, i2d_root, outdir, ["F090W"])) == 1


def test_run_verify_require_does_not_narrow_the_manifests(tmp_path, monkeypatch):
    """The manifests stay a whole-program inventory however the scope is set.

    Scoping the exit code must not drop the other filters' rows: the manifest
    is disk truth for every filter, and losing rows is how the earlier
    truncated-file bug stayed invisible.
    """
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02103_00001_nrca1_cal.fits", 1000)
    outdir = tmp_path / "out"
    _patch_mast(monkeypatch, SCOPE_ROWS)

    download.run_verify(_scope_args(cal_root, i2d_root, outdir, ["F090W"]))

    with open(outdir / "stage2_cal_manifest.csv", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert {r["filter_name"] for r in rows} == {"F090W", "F200W"}
    assert len(rows) == len(SCOPE_ROWS)


SCOPE_KIND_ROWS = [
    # F090W cal: present, exact size
    ("jw02731001001_02103_00001_nrca1_cal.fits", "CAL", "F090W", "mast:JWST/product/a1", 1000, "obs-f090w"),
    # F090W per-exposure i2d: deliberately not downloaded
    ("jw02731001001_02103_00001_nrca1_i2d.fits", "I2D", "F090W", "mast:JWST/product/a2", 500, "obs-f090w"),
    # F200W cal: present, exact size
    ("jw02731001001_02105_00001_nrca1_cal.fits", "CAL", "F200W", "mast:JWST/product/b1", 2000, "obs-f200w"),
]


def _scope_kind_fixture(tmp_path):
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02103_00001_nrca1_cal.fits", 1000)
    _place(cal_root, "jw02731001001_02105_00001_nrca1_cal.fits", 2000)
    return download.verify_products(
        _products_table(SCOPE_KIND_ROWS), cal_root, i2d_root,
        filter_name=["F090W", "F200W"],
    )


def test_require_kind_relaxes_a_partially_downloaded_filter(tmp_path):
    """F090W has all its cal files but only the combined i2d, so gating on
    F090W alone can never pass; scoping the kind can."""
    _cal, i2d = _scope_kind_fixture(tmp_path)
    gating, relaxed = i2d.split_by_scope({"F090W"}, {"CAL"})
    assert relaxed and all(r.item.kind == "I2D" for r in relaxed)
    assert not [r for r in gating if r.status == "missing"]


def test_require_kind_is_orthogonal_to_the_filter_scope(tmp_path):
    _cal, i2d = _scope_kind_fixture(tmp_path)
    # kind scope alone relaxes the F090W i2d but keeps the filter scope open
    gating, relaxed = i2d.split_by_scope(None, {"CAL"})
    assert any(r.item.filter_name == "F090W" for r in relaxed)
    # filter scope alone is not enough
    gating2, _ = i2d.split_by_scope({"F090W"}, None)
    assert any(r.status == "missing" for r in gating2)


def test_require_kind_absent_still_relaxes_nothing_when_kind_matches(tmp_path):
    _cal, i2d = _scope_kind_fixture(tmp_path)
    gating, _relaxed = i2d.split_by_scope({"F090W"}, {"I2D"})
    assert any(r.status == "missing" for r in gating)


def test_run_verify_require_kind_exits_zero_end_to_end(tmp_path, monkeypatch):
    cal_root = tmp_path / "cal"
    i2d_root = tmp_path / "i2d"
    _place(cal_root, "jw02731001001_02103_00001_nrca1_cal.fits", 1000)
    _place(cal_root, "jw02731001001_02105_00001_nrca1_cal.fits", 2000)
    outdir = tmp_path / "out"
    _patch_mast(monkeypatch, SCOPE_KIND_ROWS)
    argv = [
        "--verify",
        "--input-dir", str(cal_root),
        "--i2d-dir", str(i2d_root),
        "--outdir", str(outdir),
        "--require", "F090W",
        "--require-kind", "CAL",
    ]
    assert download.run_verify(download.parse_args(argv)) == 0
    # dropping the kind scope brings the absent F090W i2d back into the gate
    assert download.run_verify(_scope_args(cal_root, i2d_root, outdir, ["F090W"])) == 1


def test_cli_verify_parses_require_flags():
    from jwst_stack import cli

    args = cli.build_parser().parse_args(
        ["verify", "--require", "F090W", "--require-kind", "CAL"]
    )
    assert args.require == ["F090W"]
    assert args.require_kind == ["CAL"]


def test_cli_download_parses_stage2_flags():
    """The cli.py download parser is separate from download.parse_args.

    A flag added to only one of the two is invisible to the other, so the
    Stage-2 filter/kind selection has to be asserted on *both* parsers.
    """
    from jwst_stack import cli

    args = cli.build_parser().parse_args(
        ["download", "--stage2", "--filters", "F187N", "--stage2-kinds", "cal,combined_i2d"]
    )
    assert args.stage2 is True
    assert args.filters == ["F187N"]
    assert args.stage2_kinds == "cal,combined_i2d"
    # the default must keep the historical F200W, all-three-categories plan
    default = cli.build_parser().parse_args(["download", "--stage2"])
    assert default.filters == ""
    assert default.stage2_kinds == "cal,combined_i2d,subset_i2d"
    assert download.resolve_stage2_filters(default) == ["F200W"]
    assert download._parse_stage2_kinds(default.stage2_kinds) == [
        "cal",
        "combined_i2d",
        "subset_i2d",
    ]


def test_run_stage2_tolerates_missing_stage2_kinds_attr():
    """A Namespace without the attribute must not explode.

    Guards the two-parser split: run_stage2 is reachable from cli.py, whose
    parser is not download.parse_args.
    """
    prods = _two_filter_stage2_table()
    args = download.argparse.Namespace(
        proposal_id=2731,
        input_dir="/tmp/sk-data",
        i2d_dir="/tmp/sk-i2d",
        outdir="/tmp/sk-out",
        manifest=None,
        filters=["F187N"],
        plan_only=True,
        yes=True,
    )
    original = {n: getattr(download, n) for n in
                ("query_obs_table", "query_products", "check_disk", "_load_grid_shape")}
    try:
        download.query_obs_table = lambda **kw: _obs_table()
        download.query_products = lambda rows: prods
        download.check_disk = lambda *a, **k: (1e12, True, "ok")
        download._load_grid_shape = lambda outdir: (100, 100)
        assert download.run_stage2(args) == 0
    finally:
        for name, fn in original.items():
            setattr(download, name, fn)
