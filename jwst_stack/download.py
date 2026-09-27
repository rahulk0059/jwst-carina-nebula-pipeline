"""MAST product discovery and size-verified, resumable downloads.

The MAST query functions (``query_obs_table``, ``query_products``) are
module-level and monkeypatchable so all unit tests stay off the network.
Downloads are skipped when a file exists with the exact expected size, and
resumed with an HTTP ``Range`` request when a partial file is present.

``--yes`` is always optional; without it the CLI asks for interactive
confirmation before any bytes are transferred.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
import time
import urllib.request
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
from astropy.table import Table

PROPOSAL_ID = 2731
OBS_COLLECTION = "JWST"
DEFAULT_DATA_ROOT = r"C:\data\jwst_cal\mastDownload\JWST"
DEFAULT_I2D_ROOT = r"C:\data\jwst_i2d\mastDownload\JWST"

DETECTOR_RE = re.compile(
    r"_(?P<detector>nrc[ab][1-4]|nrcalong|nrcblong)_(?P<kind>cal|i2d)\.fits$",
    re.IGNORECASE,
)
_VISIT_RE = re.compile(r"^jw\d{5}\d{3}(?P<visit>\d{3})_")
_EXPOSURE_RE = re.compile(r"_(\d+)_(?P<detector>nrc[ab][1-4]|nrcblong|nrcalong)_")

MANIFEST_FIELDS = [
    "obs_id",
    "filename",
    "kind",
    "filter_name",
    "pupil",
    "detector",
    "visit",
    "exposure",
    "size",
    "status",
    "local_path",
]


@dataclass
class DownloadItem:
    """One file planned or performed for download."""

    filename: str
    url: str
    expected_size: int
    filter_name: str
    pupil: str
    detector: str
    visit: int
    exposure: int
    kind: str
    obs_id: str
    status: str = "pending"
    local_path: str = ""


def parse_product_name(filename: str) -> tuple[str, str] | None:
    """Return ``(detector, kind)`` for a JWST cal/i2d file, else None."""
    match = DETECTOR_RE.search(filename)
    if not match:
        return None
    return match.group("detector").lower(), match.group("kind").lower()


def split_filter_pupil(filters: str) -> tuple[str, str]:
    """Split a MAST ``filters`` string like ``F444W;F470N`` into (filter, pupil)."""
    parts = str(filters).split(";")
    f = parts[0].strip().upper()
    p = parts[1].strip().upper() if len(parts) > 1 else "CLEAR"
    return f, p


def visit_and_exposure_from_name(filename: str) -> tuple[int, int]:
    """Best-effort (visit, exposure) from the file/obs name structure."""
    vm = _VISIT_RE.match(filename)
    visit = int(vm.group("visit")) if vm else -1
    em = _EXPOSURE_RE.search(filename)
    exposure = int(em.group(1)) if em else -1
    return visit, exposure


def query_obs_table(
    proposal_id: int = PROPOSAL_ID,
    obs_collection: str = OBS_COLLECTION,
    dataproduct_type: str = "image",
) -> Table:
    """Query the MAST observation table (module-level, injectable for tests)."""
    from astroquery.mast import Observations

    return Observations.query_criteria(
        proposal_id=proposal_id,
        obs_collection=obs_collection,
        dataproduct_type=dataproduct_type,
    )


def query_products(obs_rows: Table) -> Table:
    """Query the MAST product list for the given observation rows."""
    from astroquery.mast import Observations

    return Observations.get_product_list(obs_rows)


MAST_SCHEME = "https://mast.stsci.edu/api/v0.1/Download/file?uri="


def mast_download_url(data_uri: str) -> str:
    """Public HTTPS URL for a MAST ``dataURI``."""
    return f"{MAST_SCHEME}{data_uri}"


def is_nircam_obs(obs_id: str) -> bool:
    return "nircam" in obs_id.lower()


def _make_item(row, short: str) -> DownloadItem:
    """Turn one MAST product row into a :class:`DownloadItem`.

    ``short`` is the uppercase subclass (``CAL``, ``I2D``, ...).  The
    combined i2d mosaic has no detector/exposure components; its detector is
    labelled ``combined``.
    """
    filename = str(row["productFilename"])
    parsed = parse_product_name(filename)
    if parsed is None:
        detector = "combined"
        visit = exposure = 0
    else:
        detector, _ = parsed
        visit, exposure = visit_and_exposure_from_name(filename)
    filter_name, pupil = split_filter_pupil(row["filters"])
    uri = row["dataURI"]
    url = mast_download_url(uri) if uri else ""
    size = int(row["size"]) if row["size"] is not None else 0
    return DownloadItem(
        filename=filename,
        url=url,
        expected_size=size,
        filter_name=filter_name,
        pupil=pupil,
        detector=detector,
        visit=visit,
        exposure=exposure,
        kind=short,
        obs_id=str(row["obs_id"]),
    )


def build_items(
    products: Table,
    want_kinds: set[str],
    data_root: str | Path,
    filter_names: set[str] | None = None,
    detectors: set[str] | None = None,
) -> list[DownloadItem]:
    """Project a MAST product table onto :class:`DownloadItem` rows."""
    data_root = Path(data_root)
    items: list[DownloadItem] = []
    for row in products:
        filename = str(row["productFilename"])
        kind = "" if row["productSubGroupDescription"] is np.ma.masked else str(row["productSubGroupDescription"])
        kind = kind.upper()
        if kind not in want_kinds:
            continue
        parsed = parse_product_name(filename)
        if parsed is None:
            continue
        detector, _ = parsed
        if detectors is not None and detector not in detectors:
            continue
        filter_name, _ = split_filter_pupil(row["filters"])
        if filter_names is not None and filter_name not in filter_names:
            continue
        items.append(_make_item(row, kind))
    items.sort(key=lambda i: (i.filter_name, i.visit, i.detector, i.exposure, i.kind))
    return items


def summarize_items(items: list[DownloadItem]) -> str:
    lines = [f"{'filter':>8} {'kind':>4} {'visit':>5} {'detector':>9} {'files':>6} {'total MB':>12}"]
    groups: dict[tuple, list[DownloadItem]] = {}
    for it in items:
        groups.setdefault((it.filter_name, it.kind, it.visit, it.detector), []).append(it)
    for key in sorted(groups):
        f, kind, visit, det = key
        group = groups[key]
        lines.append(
            f"{f:>8} {kind:>4} {visit:>5} {det:>9} {len(group):>6} "
            f"{sum(i.expected_size for i in group) / 1e6:>12.1f}"
        )
    total = sum(i.expected_size for i in items)
    lines.append(f"\nTotal: {len(items)} files, {total/1e6:.1f} MB")
    return "\n".join(lines)


def is_combined_i2d(filename: str) -> bool:
    """True for the single all-visit depth mosaic (no exposure/detector ids)."""
    return filename.lower().endswith("_i2d.fits") and parse_product_name(filename) is None


def build_i2d_items(
    products: Table,
    data_root: str | Path,
    filter_name: str | Sequence[str] = "F200W",
) -> list[DownloadItem]:
    """Every i2d product, *including* the combined all-visit mosaic.

    ``build_items`` cannot be used for i2d because it drops products whose
    name has no detector/exposure components -- which is exactly the combined
    mosaic, the one file whose size is easiest to get wrong.

    ``filter_name`` may be one filter or several.
    """
    if isinstance(filter_name, str):
        wanted = {filter_name.upper()}
    else:
        wanted = {str(f).upper() for f in filter_name if str(f).strip()}
    items: list[DownloadItem] = []
    for row in products:
        filename = str(row["productFilename"])
        kind = (
            ""
            if row["productSubGroupDescription"] is np.ma.masked
            else str(row["productSubGroupDescription"])
        )
        if kind.upper() != "I2D":
            continue
        fname, _ = split_filter_pupil(row["filters"])
        if fname.upper() not in wanted:
            continue
        items.append(_make_item(row, "I2D"))
    items.sort(key=lambda i: (is_combined_i2d(i.filename), i.detector, i.visit, i.exposure))
    return items


def plan_stage2(
    products: Table,
    data_root: str | Path,
    i2d_root: str | Path,
    filter_names: str | Sequence[str] = "F200W",
) -> dict:
    """Classify one filter's products into the Stage-2 selection and disk status.

    Returns a dict with the planned categories: every cal file, the combined
    i2d mosaic, and the per-exposure i2d validation subset (nrca1 for all
    visits plus all eight SW detectors for visit 1).  Files already present
    with a matching size are flagged ``existing``.

    ``filter_names`` may be one filter or several, and defaults to F200W only
    when the caller does not say.  It must be honoured: the categories, the
    reported sizes and the printed labels are all filter-specific, so silently
    substituting a default here mislabels the whole plan.  The resolved filters
    are returned as ``plan["filters"]`` for the formatter to label with.
    """
    if isinstance(filter_names, str):
        wanted = [filter_names.upper()]
    else:
        wanted = [str(f).upper() for f in filter_names if str(f).strip()]
    if not wanted:
        raise ValueError("plan_stage2 needs at least one filter name")

    data_root = Path(data_root)
    i2d_root = Path(i2d_root)
    all_cal = build_items(products, {"CAL"}, data_root, filter_names=set(wanted))

    all_i2d = build_i2d_items(products, i2d_root, filter_name=wanted)
    combined = [i for i in all_i2d if is_combined_i2d(i.filename)]
    perexposure = [i for i in all_i2d if not is_combined_i2d(i.filename)]

    subset = [
        it
        for it in perexposure
        if (it.detector == "nrca1" and it.visit in {1, 2, 3, 4}) or (it.visit == 1)
    ]

    cal_index = _index_existing(data_root)
    i2d_index = _index_existing(i2d_root)

    def status(item: DownloadItem, index: dict[str, Path]) -> str:
        if _is_existing(index, item):
            item.status = "existing"
            return "existing"
        item.status = "new"
        return "new"

    categories = {
        "cal": [(it, status(it, cal_index)) for it in all_cal],
        "combined_i2d": [(it, status(it, i2d_index)) for it in combined],
        "subset_i2d": [(it, status(it, i2d_index)) for it in subset],
    }
    return {
        "categories": categories,
        "filters": wanted,
        "all_cal": all_cal,
        "all_i2d": all_i2d,
        "combined": combined,
        "subset": subset,
    }


_STAGE2_CATEGORIES = (
    ("cal", "cal files (all)"),
    ("combined_i2d", "combined i2d mosaic"),
    ("subset_i2d", "per-exposure i2d validation subset"),
)


def _parse_stage2_kinds(spec: str) -> list[str]:
    """Validate and order a ``--stage2-kinds`` selection."""
    valid = {k for k, _ in _STAGE2_CATEGORIES}
    wanted = [k.strip() for k in str(spec).split(",") if k.strip()]
    bad = [k for k in wanted if k not in valid]
    if bad:
        raise ValueError(
            f"unknown --stage2-kinds {bad}; valid: {sorted(valid)}"
        )
    if not wanted:
        raise ValueError("--stage2-kinds must name at least one category")
    return [k for k, _ in _STAGE2_CATEGORIES if k in wanted]


def format_stage2_plan(
    plan: dict,
    scratch_estimate_gb: float | None = None,
    filters: Sequence[str] | None = None,
    kinds: Sequence[str] | None = None,
) -> str:
    if filters is None:
        filters = plan.get("filters") or ["F200W"]
    label_filter = ",".join(str(f).upper() for f in filters)
    if kinds is None:
        kinds = [k for k, _ in _STAGE2_CATEGORIES]
    kinds = [k for k in kinds if k in dict(_STAGE2_CATEGORIES)]

    lines = []
    grand_existing = grand_pending = 0
    for kind, desc in _STAGE2_CATEGORIES:
        if kind not in kinds:
            continue
        rows = plan["categories"][kind]
        existing = [it for it, s in rows if s == "existing"]
        pending = [it for it, s in rows if s != "existing"]
        grand_existing += len(existing)
        grand_pending += len(pending)
        lines.append(
            f"{f'{label_filter} {desc}':<44} {len(rows):>4} planned  "
            f"{len(existing):>4} on disk/verified  {len(pending):>4} to download  "
            f"{sum(i.expected_size for i in pending) / 1e9:6.2f} GB"
        )
    total_pending = sum(
        it.expected_size
        for kind in kinds
        for it, s in plan["categories"][kind]
        if s != "existing"
    )
    lines.append(f"\n  total to download : {grand_pending} files, {total_pending/1e9:.2f} GB")
    if scratch_estimate_gb is not None:
        lines.append(f"  est. mosaic scratch: {scratch_estimate_gb:.1f} GB peak (Stage 4)")
    return "\n".join(lines)


def estimate_mosaic_scratch_gb(
    naxis: tuple[int, int],
    n_frames: int,
    overlay_naxis: int = 2050,
    bytes_per_px_output: int = 10,
    bytes_per_px_overlay: int = 10,
) -> float:
    """Peak temporary disk (GB) for a tiled mosaic build.

    ``bytes_per_px_output`` covers the mosaic SCI (float32) + WEIGHT (float32)
    + COVERAGE (uint16) arrays; ``bytes_per_px_overlay`` covers one reprojected
    footprint overlay (SCI + ERR float32 + coverage uint16) while the final
    arrays exist concurrently.
    """
    n_px = naxis[0] * naxis[1]
    overlays = n_frames * overlay_naxis * overlay_naxis
    return (n_px * bytes_per_px_output + overlays * bytes_per_px_overlay) / 1e9


def check_disk(
    dest_root: str | Path, total_bytes: int, safety: float = 1.1
) -> tuple[float, bool, str]:
    """Free bytes on the drive holding *dest_root* and whether it is enough."""
    dest_root = Path(dest_root)
    dest_root.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(dest_root).free
    needed = int(total_bytes * safety)
    if free >= needed:
        return float(free), True, f"{free/1e9:.1f} GB free >= {needed/1e9:.1f} GB needed"
    return (
        float(free),
        False,
        f"{free/1e9:.1f} GB free < {needed/1e9:.1f} GB needed "
        f"({total_bytes/1e6:.0f} MB planned, {safety:.0%} safety)",
    )


def _group_dir_name(filename: str) -> str:
    """MAST-ish subfolder for a product (group id): name minus kind suffix."""
    lowered = filename.lower()
    for suffix in ("_cal.fits", "_i2d.fits"):
        if lowered.endswith(suffix):
            return filename[: -len(suffix)]
    return ""


def _destination_path(root: Path, item: DownloadItem) -> Path:
    return root / _group_dir_name(item.filename) / item.filename


def _index_existing(root: Path) -> dict[str, Path]:
    """Map every file name under *root* to its path (recursive)."""
    if not root.exists():
        return {}
    return {
        p.name: p
        for p in root.rglob("*")
        if p.is_file() and p.name.lower().endswith((".fits", ".fit"))
    }


def _is_existing(index: dict[str, Path], item: DownloadItem) -> bool:
    path = index.get(item.filename)
    if path is None:
        return False
    return path.stat().st_size == item.expected_size


def file_size_matches(path: Path, expected: int) -> bool:
    if not path.exists():
        return False
    return path.stat().st_size == expected


def download_url(
    url: str,
    dest: Path,
    expected_size: int,
    resume: bool = True,
    chunk: int = 1 << 16,
    progress: Callable[[int, int], None] | None = None,
    max_seconds: float | None = None,
    report_every: int = 1 << 28,
) -> str:
    """Download *url* to *dest*, skipping/continuing partial files.

    Returns one of ``skipped``, ``ok``, ``resumed``, ``size_mismatch``,
    ``error``.

    ``progress`` is called as ``(bytes_written, expected_size)`` every
    ``report_every`` bytes.  It exists because a large single product (the
    combined i2d is 5.4 GB against 117 MB for a cal file) otherwise emits no
    output whatsoever for the whole transfer, which makes a slow transfer
    indistinguishable from a hang.  ``max_seconds`` is a wall-clock cap on the
    transfer, checked between chunks, so a pathological trickle is abandoned
    and reported instead of running forever.
    """
    if file_size_matches(dest, expected_size):
        return "skipped"
    if dest.exists():
        current = dest.stat().st_size
        if current > expected_size:
            dest.unlink()
            current = 0
        if current == expected_size:
            return "skipped"
    else:
        current = 0

    headers = {}
    if resume and current > 0:
        headers["Range"] = f"bytes={current}-"

    started = time.monotonic()
    next_report = current
    try:
        request = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(request, timeout=120) as response:
            dest.parent.mkdir(parents=True, exist_ok=True)
            # A server that answers 200 to a Range request is sending the whole
            # body, so the partial must be truncated.  It cannot be done with
            # seek(0) on a handle opened in "ab": O_APPEND forces every write
            # to the end regardless of the file position, which would append
            # the full body onto the partial and leave an oversized corrupt
            # file behind.
            restarting = response.status == 200 and current > 0 and resume
            if restarting:
                current = 0
            mode = "ab" if (resume and current > 0) else "wb"
            with open(dest, mode) as handle:
                while True:
                    if max_seconds is not None and time.monotonic() - started > max_seconds:
                        return (
                            f"error:abandoned after {time.monotonic() - started:.0f}s "
                            f"at {dest.stat().st_size}/{expected_size} bytes"
                        )
                    block = response.read(chunk)
                    if not block:
                        break
                    handle.write(block)
                    written = handle.tell()
                    if progress is not None and written - next_report >= report_every:
                        next_report = written
                        progress(written, expected_size)
    except Exception as exc:  # pragma: no cover - network resilience
        return f"error:{exc}"

    if dest.stat().st_size == expected_size:
        return "resumed" if current > 0 and "Range" in headers else "ok"
    return "size_mismatch"


def item_failed(item: DownloadItem) -> bool:
    """True when *item* did not end up complete and size-verified on disk.

    Only ``ok`` (and ``skipped``, for an item with no URL) counts as success.
    Anything else -- a network ``error:``, a ``size_mismatch``, or an
    interrupted transfer -- must gate, otherwise a partial download reports
    itself complete.
    """
    return item.status not in {"ok", "skipped"}


def _print_line(message: str) -> None:
    """Print and flush.

    Flushing matters: when stdout is redirected to a log file it is block
    buffered, so an unflushed progress line for a 5 GB transfer would not
    appear until the process exited -- which is exactly the "no output at
    all" symptom that made a slow download look like a hang.
    """
    print(message, flush=True)


def download_items(
    items: list[DownloadItem],
    dest_root: str | Path,
    manifest_path: str | Path,
    progress: Callable[[str], None] | None = None,
) -> list[DownloadItem]:
    """Download *items* into *dest_root* and record status in the manifest."""
    emit = progress or _print_line
    dest_root = Path(dest_root)
    dest_root.mkdir(parents=True, exist_ok=True)
    for item in items:
        dest = _destination_path(dest_root, item)
        dest.parent.mkdir(parents=True, exist_ok=True)
        item.local_path = str(dest)
        if item.url:
            def report(written: int, total: int, _name=item.filename) -> None:
                pct = 100.0 * written / total if total else 0.0
                emit(
                    f"  {_name:45s} {written/1e6:9.1f}/{total/1e6:.1f} MB "
                    f"({pct:5.1f}%)"
                )

            status = download_url(
                item.url, dest, item.expected_size, resume=True, progress=report
            )
        else:
            status = "skipped"
        item.status = status
        if status.startswith("ok") or status.startswith("resumed"):
            item.status = "ok"
        emit(f"{item.filename:45s} {item.expected_size/1e6:8.1f} MB  {item.status}")
    write_manifest(items, manifest_path)
    return items


def write_manifest(items: list[DownloadItem], manifest_path: str | Path) -> Path:
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        for item in items:
            writer.writerow(
                {
                    "obs_id": item.obs_id,
                    "filename": item.filename,
                    "kind": item.kind,
                    "filter_name": item.filter_name,
                    "pupil": item.pupil,
                    "detector": item.detector,
                    "visit": item.visit,
                    "exposure": item.exposure,
                    "size": item.expected_size,
                    "status": item.status,
                    "local_path": item.local_path,
                }
            )
    return manifest_path


def read_manifest(manifest_path: str | Path) -> list[dict]:
    manifest_path = Path(manifest_path)
    if not manifest_path.exists():
        return []
    with open(manifest_path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="jwst_stack download")
    parser.add_argument("--input-dir", default=DEFAULT_DATA_ROOT, help="destination root")
    parser.add_argument("--i2d-dir", default=DEFAULT_I2D_ROOT, help="i2d destination root")
    parser.add_argument("--outdir", default="out")
    parser.add_argument("--manifest", default=None, help="manifest CSV path")
    parser.add_argument("--filters", default="", nargs="+", help="e.g. F200W")
    parser.add_argument("--detectors", default="", nargs="+", help="e.g. nrca1 nrcblong")
    parser.add_argument(
        "--products",
        default="cal,i2d",
        help="comma-separated: cal,i2d",
    )
    parser.add_argument("--proposal-id", type=int, default=PROPOSAL_ID)
    parser.add_argument(
        "--stage2",
        action="store_true",
        help=(
            "pre-selected Stage-2 set: all cal files, the combined i2d mosaic "
            "and the per-exposure i2d validation subset, for each --filters "
            "filter (default: F200W)"
        ),
    )
    parser.add_argument(
        "--stage2-kinds",
        default="cal,combined_i2d,subset_i2d",
        help=(
            "comma-separated subset of the Stage-2 categories to plan and "
            "download: cal, combined_i2d, subset_i2d. The per-exposure i2d "
            "subset_i2d is only needed for per-exposure validation; Stage-4 "
            "mosaic validation uses combined_i2d alone."
        ),
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="compare disk against the MAST plan and rebuild the manifests",
    )
    parser.add_argument(
        "--filter",
        dest="filter_name",
        nargs="+",
        default=None,
        help=(
            "filters to verify (default: every NIRCam filter in the MAST plan). "
            "Must cover every filter present under the data root -- a file whose "
            "filter is not in the plan is reported as 'extra'."
        ),
    )
    parser.add_argument(
        "--require",
        nargs="+",
        default=None,
        metavar="FILTER",
        help=(
            "only these filters gate the exit code (e.g. --require F090W). "
            "Absent products for other filters are still listed in the report "
            "but do not fail the run. The scan and the manifests always cover "
            "every filter, so scoping never narrows the disk-truth inventory."
        ),
    )
    parser.add_argument(
        "--require-kind",
        nargs="+",
        default=None,
        metavar="KIND",
        help=(
            "only these product kinds gate the exit code (CAL, I2D). Combine "
            "with --require when a filter was only partially downloaded, e.g. "
            "--require F090W --require-kind CAL when the per-exposure i2d "
            "products were deliberately skipped"
        ),
    )
    parser.add_argument(
        "--no-rebuild",
        action="store_true",
        help="report only, do not rewrite the manifests",
    )
    parser.add_argument("--plan-only", action="store_true", help="print the plan, do not download")
    parser.add_argument("--yes", action="store_true", help="skip the interactive prompt")
    return parser.parse_args(argv)


def resolve_stage2_filters(args: argparse.Namespace) -> list[str]:
    """Filters the Stage-2 plan covers, upper-cased.

    ``--filters`` is honoured; F200W is the fallback only when no filter was
    requested, never an override of one that was.
    """
    requested = getattr(args, "filters", None) or []
    wanted = [str(f).strip().upper() for f in requested if str(f).strip()]
    return wanted or ["F200W"]


def run_stage2(args: argparse.Namespace, confirm: Callable[[str], bool] | None = None) -> int:
    obs_rows = query_obs_table(proposal_id=args.proposal_id)
    nircam = obs_rows[[is_nircam_obs(o) for o in obs_rows["obs_id"]]]
    products = query_products(nircam)
    filters = resolve_stage2_filters(args)
    kinds = _parse_stage2_kinds(
        getattr(args, "stage2_kinds", None) or "cal,combined_i2d,subset_i2d"
    )
    plan = plan_stage2(products, args.input_dir, args.i2d_dir, filters)

    grid_shape = _load_grid_shape(Path(args.outdir))
    scratch = estimate_mosaic_scratch_gb(grid_shape, len(plan["all_cal"]))

    print(f"\nStage-2 download plan ({', '.join(filters)}):")
    print(format_stage2_plan(plan, scratch, filters, kinds))

    pending = [
        it
        for kind in kinds
        for it, status in plan["categories"][kind]
        if status != "existing"
    ]
    total_bytes = sum(i.expected_size for i in pending)
    free, ok, message = check_disk(args.input_dir, total_bytes)
    print(f"\nDisk ({args.input_dir}): {message}")
    if not ok:
        return 1
    if args.plan_only:
        print("plan-only: nothing downloaded")
        return 0

    if not args.yes:
        if confirm is None:
            reply = input("Proceed with download? [y/N] ")
        else:
            reply = str(confirm(message)).strip()
        if not (reply[:1].lower() == "y"):
            print("cancelled")
            return 0

    cal_pending = (
        [it for it, s in plan["categories"]["cal"] if s != "existing"]
        if "cal" in kinds
        else []
    )
    i2d_pending = [
        it
        for kind in ("combined_i2d", "subset_i2d")
        if kind in kinds
        for it, s in plan["categories"][kind]
        if s != "existing"
    ]
    downloaded: list[DownloadItem] = []
    if cal_pending:
        print(f"\ncal files -> {args.input_dir}")
        downloaded += download_items(
            cal_pending,
            args.input_dir,
            Path(args.manifest) if args.manifest else Path(args.outdir) / "stage2_cal_manifest.csv",
        )
    if i2d_pending:
        print(f"\ni2d files -> {args.i2d_dir}")
        downloaded += download_items(
            i2d_pending,
            args.i2d_dir,
            Path(args.manifest) if args.manifest else Path(args.outdir) / "stage2_i2d_manifest.csv",
        )
    failed = [it for it in downloaded if item_failed(it)]
    if failed:
        print(
            f"\nStage-2 download INCOMPLETE: {len(failed)} of {len(downloaded)} "
            f"file(s) failed. Re-run the same command to resume."
        )
        for it in failed[:10]:
            print(f"  {it.filename:45s} {it.status}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more")
        return 1
    print("\nStage-2 download complete")
    return 0


def _load_grid_shape(outdir: str | Path) -> tuple[int, int]:
    """Grid shape from out/grid.fits, falling back to the known F200W value."""
    grid_path = Path(outdir) / "grid.fits"
    if grid_path.exists():
        from jwst_stack import grid as grid_mod

        _, shape = grid_mod.load_grid(grid_path)
        return shape
    return (15895, 22130)


def run_download(args: argparse.Namespace, confirm: Callable[[str], bool] | None = None) -> int:
    want_kinds = {k.upper() for k in str(args.products).split(",")}
    obs_rows = query_obs_table(proposal_id=args.proposal_id)
    nircam = obs_rows[[is_nircam_obs(o) for o in obs_rows["obs_id"]]]
    products = query_products(nircam)
    filter_names = {f.upper() for f in args.filters} if args.filters else None
    detectors = {d.lower() for d in args.detectors} if args.detectors else None
    items = build_items(
        products, want_kinds, args.input_dir, filter_names, detectors
    )
    if not items:
        print("no products matched the requested filters/detectors")
        return 1

    print("\nPlanned downloads:")
    print(summarize_items(items))

    total_bytes = sum(i.expected_size for i in items)
    free, ok, message = check_disk(args.input_dir, total_bytes)
    print(f"Disk ({args.input_dir}): {message}")
    if not ok:
        return 1

    if args.plan_only:
        print("plan-only: nothing downloaded")
        return 0

    if not args.yes:
        if confirm is None:
            reply = input("Proceed with download? [y/N] ")
        else:
            reply = str(confirm(message)).strip()
        if not (reply[:1].lower() == "y"):
            print("cancelled")
            return 0

    manifest = Path(args.manifest) if args.manifest else Path(args.outdir) / "download_manifest.csv"
    downloaded = download_items(items, args.input_dir, manifest)
    print(f"manifest: {manifest}")
    failed = [it for it in downloaded if item_failed(it)]
    if failed:
        print(
            f"download INCOMPLETE: {len(failed)} of {len(downloaded)} file(s) "
            f"failed. Re-run the same command to resume."
        )
        for it in failed[:10]:
            print(f"  {it.filename:45s} {it.status}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more")
        return 1
    return 0


@dataclass
class VerifyRow:
    """One expected product compared against what is actually on disk."""

    item: DownloadItem
    status: str
    actual_size: int
    local_path: str = ""


#: Statuses that gate the exit code whatever the filter scope, because they mean
#: a file on disk is corrupt or unexpected rather than merely not downloaded.
_ALWAYS_GATING = frozenset({"size_mismatch", "extra"})


@dataclass
class VerifyReport:
    """Disk-truth result for one destination root."""

    label: str
    root: Path
    rows: list[VerifyRow] = field(default_factory=list)

    def count(self, status: str) -> int:
        return sum(1 for r in self.rows if r.status == status)

    @property
    def problems(self) -> list[VerifyRow]:
        return [r for r in self.rows if r.status != "ok"]

    def split_by_scope(
        self,
        scope: set[str] | None,
        kinds: set[str] | None = None,
    ) -> tuple[list[VerifyRow], list[VerifyRow]]:
        """Return ``(gating, out_of_scope)`` problems for an optional scope.

        A row is in scope when it matches *scope* (filters) and *kinds* (product
        kinds such as ``CAL`` / ``I2D``); either may be ``None`` to leave that
        dimension unrestricted.  With no scope at all every problem gates, which
        is the whole-program behaviour.

        Absent products for an out-of-scope filter or kind do not gate, because
        a product that was never downloaded is not a disk problem.  This is what
        lets a partially-downloaded filter exit clean, e.g. F090W has all 160
        cal files but only the combined i2d, so ``--require F090W
        --require-kind CAL`` is satisfied.

        ``size_mismatch`` and ``extra`` always gate regardless of scope: a
        truncated or unexpected file under a data root is a real anomaly
        whichever filter or kind is being checked.
        """
        if scope is None and kinds is None:
            return list(self.problems), []
        gating: list[VerifyRow] = []
        relaxed: list[VerifyRow] = []
        for row in self.problems:
            if row.status in _ALWAYS_GATING:
                gating.append(row)
                continue
            in_filters = scope is None or str(row.item.filter_name or "").upper() in scope
            in_kinds = kinds is None or str(row.item.kind or "").upper() in kinds
            (gating if in_filters and in_kinds else relaxed).append(row)
        return gating, relaxed

    @property
    def total_bytes(self) -> int:
        return sum(r.actual_size for r in self.rows)

    def expected_bytes(self) -> int:
        return sum(r.item.expected_size for r in self.rows)


def _classify(
    items: list[DownloadItem], index: dict[str, Path], root: Path
) -> list[VerifyRow]:
    """Compare planned *items* against the files present under *root*."""
    rows: list[VerifyRow] = []
    seen: set[str] = set()
    for item in items:
        path = index.get(item.filename)
        if path is None:
            status, actual = "missing", 0
            local = ""
        else:
            seen.add(item.filename)
            actual = path.stat().st_size
            local = str(path)
            if item.expected_size and actual == item.expected_size:
                status = "ok"
            else:
                status = "size_mismatch"
        row = VerifyRow(item=item, status=status, actual_size=actual, local_path=local)
        row.item.status = status
        row.item.local_path = local
        rows.append(row)

    for filename, path in sorted(index.items()):
        if filename in seen:
            continue
        extra = DownloadItem(
            filename=filename,
            url="",
            expected_size=0,
            filter_name="",
            pupil="",
            detector="",
            visit=0,
            exposure=0,
            kind="",
            obs_id="",
            status="extra",
            local_path=str(path),
        )
        rows.append(
            VerifyRow(
                item=extra,
                status="extra",
                actual_size=path.stat().st_size,
                local_path=str(path),
            )
        )
    rows.sort(key=lambda r: (r.status != "ok", r.status, r.item.filename))
    return rows


def nircam_filters(products: Table) -> list[str]:
    """Every distinct NIRCam filter named in a MAST product list, sorted."""
    names: set[str] = set()
    for row in products:
        fname, _ = split_filter_pupil(row["filters"])
        if fname:
            names.add(fname.upper())
    return sorted(names)


def verify_products(
    products: Table,
    data_root: str | Path,
    i2d_root: str | Path,
    filter_name: str | Sequence[str] = "F200W",
) -> tuple[VerifyReport, VerifyReport]:
    """Compare every planned cal / i2d product with the files on disk.

    The manifests record what was *attempted*; this reports what is actually
    there, so a stale manifest can never hide a truncated or missing file.

    ``filter_name`` may be a single filter or several.  It must cover every
    filter actually present under the data root: a file whose filter is absent
    from the plan is classified ``extra``, so a single-filter plan over a
    multi-filter root reports every other filter's files as problems and the
    non-zero exit can no longer gate anything.
    """
    if isinstance(filter_name, str):
        filter_names = {filter_name.upper()}
    else:
        filter_names = {str(f).upper() for f in filter_name if str(f).strip()}

    cal_items = build_items(products, {"CAL"}, data_root, filter_names=filter_names)
    i2d_items = build_i2d_items(products, i2d_root, filter_name=filter_names)

    data_root = Path(data_root)
    i2d_root = Path(i2d_root)
    cal_report = VerifyReport(
        label="cal", root=data_root, rows=_classify(cal_items, _index_existing(data_root), data_root)
    )
    i2d_report = VerifyReport(
        label="i2d", root=i2d_root, rows=_classify(i2d_items, _index_existing(i2d_root), i2d_root)
    )
    return cal_report, i2d_report


def write_verify_manifest(report: VerifyReport, manifest_path: str | Path) -> Path:
    """Rewrite a manifest so it matches the files actually on disk."""
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        for row in report.rows:
            item = row.item
            writer.writerow(
                {
                    "obs_id": item.obs_id,
                    "filename": item.filename,
                    "kind": item.kind,
                    "filter_name": item.filter_name,
                    "pupil": item.pupil,
                    "detector": item.detector,
                    "visit": item.visit,
                    "exposure": item.exposure,
                    "size": row.actual_size,
                    "status": row.status,
                    "local_path": row.local_path,
                }
            )
    return manifest_path


def _verify_problem_table(title: str, rows: list[VerifyRow]) -> list[str]:
    lines = [f"  {title}"]
    if not rows:
        lines.append("    none")
        return lines
    lines.append(f"    {'status':>14} {'expected':>13} {'actual':>13}  filename")
    for row in rows:
        expected = (
            f"{row.item.expected_size:>13,}"
            if row.item.expected_size
            else f"{'-':>13}"
        )
        lines.append(
            f"    {row.status:>14} {expected} {row.actual_size:>13,}  "
            f"{row.item.filename}"
        )
    return lines


def format_verify_report(
    reports: list[VerifyReport],
    scope: set[str] | None = None,
    kinds: set[str] | None = None,
) -> str:
    """Human-readable disk-truth summary, with every problem listed.

    *scope*/*kinds* limit which rows gate the exit code; out-of-scope absences
    are still listed, under their own heading, so scoping never hides anything.
    """
    lines: list[str] = []
    for report in reports:
        gating, relaxed = report.split_by_scope(scope, kinds)
        lines.append(f"\n{report.label}: {report.root}")
        lines.append(
            f"  {report.count('ok')} ok, {report.count('size_mismatch')} size_mismatch, "
            f"{report.count('missing')} missing, {report.count('extra')} extra "
            f"({len(report.rows)} rows, {report.total_bytes/1e9:.3f} GB on disk)"
        )
        if scope is not None or kinds is not None:
            parts = []
            if scope is not None:
                parts.append(f"filters {', '.join(sorted(scope))}")
            if kinds is not None:
                parts.append(f"kinds {', '.join(sorted(kinds))}")
            lines.append(f"  gating scope: {'; '.join(parts)}")
        lines.extend(_verify_problem_table("problems (gating)", gating))
        if relaxed:
            lines.extend(
                _verify_problem_table(
                    f"absent, outside the gating scope ({len(relaxed)} rows)", relaxed
                )
            )
    total = sum(len(r.split_by_scope(scope, kinds)[0]) for r in reports)
    relaxed_total = sum(len(r.split_by_scope(scope, kinds)[1]) for r in reports)
    lines.append("")
    if relaxed_total:
        lines.append(
            f"({relaxed_total} absent row(s) outside the gating scope are listed "
            "above but do not affect the exit code)"
        )
    if total:
        lines.append(f"VERIFY FAILED: {total} problem(s) need attention")
    else:
        lines.append("VERIFY OK: every gated file is present with the exact size")
    return "\n".join(lines)


def run_verify(args: argparse.Namespace) -> int:
    obs_rows = query_obs_table(proposal_id=args.proposal_id)
    nircam = obs_rows[[is_nircam_obs(o) for o in obs_rows["obs_id"]]]
    products = query_products(nircam)
    filter_names = getattr(args, "filter_name", None)
    if isinstance(filter_names, str):
        filter_names = [filter_names]
    if not filter_names:
        filter_names = nircam_filters(products)
        print(
            "verify: no --filter given; checking every NIRCam filter in the plan "
            f"({', '.join(filter_names)})"
        )

    # --require scopes the *exit code* only.  The scan always covers every
    # filter so the manifests stay a whole-program inventory of disk truth.
    required = getattr(args, "require", None)
    if isinstance(required, str):
        required = [required]
    scope = {str(f).upper() for f in (required or []) if str(f).strip()} or None
    kinds_arg = getattr(args, "require_kind", None)
    if isinstance(kinds_arg, str):
        kinds_arg = [kinds_arg]
    kinds = {str(k).upper() for k in (kinds_arg or []) if str(k).strip()} or None
    if scope or kinds:
        parts = []
        if scope:
            parts.append(f"filters {', '.join(sorted(scope))}")
        if kinds:
            parts.append(f"kinds {', '.join(sorted(kinds))}")
        print(
            "verify: gating on "
            f"{'; '.join(parts)}; other rows are still listed but do not "
            "affect the exit code"
        )

    cal_report, i2d_report = verify_products(
        products, args.input_dir, args.i2d_dir, filter_name=filter_names
    )
    reports = [cal_report, i2d_report]
    text = format_verify_report(reports, scope=scope, kinds=kinds)
    print(text)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    report_path = outdir / "verify_report.txt"
    report_path.write_text(text + "\n", encoding="utf-8")
    print(f"report: {report_path}")

    if not getattr(args, "no_rebuild", False):
        cal_manifest = write_verify_manifest(
            cal_report, outdir / "stage2_cal_manifest.csv"
        )
        i2d_manifest = write_verify_manifest(
            i2d_report, outdir / "stage2_i2d_manifest.csv"
        )
        print(f"rebuilt manifests: {cal_manifest}, {i2d_manifest}")

    return 1 if sum(len(r.split_by_scope(scope, kinds)[0]) for r in reports) else 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.verify:
        return run_verify(args)
    if args.stage2:
        return run_stage2(args)
    return run_download(args)


if __name__ == "__main__":
    sys.exit(main())