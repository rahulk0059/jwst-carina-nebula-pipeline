"""Command-line interface: inspect, group, stack and compare."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from jwst_stack import (
    download,
    gauge as gauge_mod,
    grid as grid_mod,
    io,
    grouping,
    plotting,
    stack as stack_mod,
)


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input-dir", default=r"C:\data\jwst_cal\mastDownload\JWST")
    parser.add_argument("--outdir", default="out")


def _add_selection(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--detector", nargs="+", default=None, help="e.g. nrca1 (repeatable)"
    )
    parser.add_argument(
        "--filter",
        dest="filter_name",
        nargs="+",
        default=None,
        help="e.g. F200W (repeatable)",
    )
    parser.add_argument(
        "--pupil",
        nargs="+",
        default=None,
        help=(
            "e.g. CLEAR or F470N (repeatable). REQUIRED whenever the selection "
            "spans more than one pupil: F444W is exposed through both CLEAR and "
            "F470N and both read FILTER=F444W, so --filter F444W alone would "
            "silently stack two bandpasses"
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="jwst_stack", description="Stack JWST NIRCam calibrated exposures."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_inspect = sub.add_parser("inspect", help="Print a per-file metadata table.")
    _add_common(p_inspect)
    _add_selection(p_inspect)

    p_group = sub.add_parser("group", help="Report footprint-overlap groups.")
    _add_common(p_group)
    _add_selection(p_group)
    p_group.add_argument("--threshold", type=float, default=0.15)
    p_group.add_argument("--sample", type=int, default=48)

    p_stack = sub.add_parser("stack", help="Align and sigma-clip stack exposures.")
    _add_common(p_stack)
    _add_selection(p_stack)
    p_stack.add_argument("--visit", default="all", help="1-4 or 'all'")
    p_stack.add_argument("--sigma", type=float, default=3.0)
    p_stack.add_argument("--iterations", type=int, default=3)
    p_stack.add_argument("--combine", choices=["median", "mean"], default="median")
    p_stack.add_argument("--method", choices=["interp", "exact"], default="interp")
    p_stack.add_argument("--pixel-scale", type=float, default=None, help="arcsec/px")
    p_stack.add_argument(
        "--no-register",
        action="store_true",
        help="skip photutils star registration (reproduce the pure-WCS stack)",
    )
    p_stack.add_argument("--fwhm-px", type=float, default=3.0)
    p_stack.add_argument(
        "--match-radius-arcsec",
        type=float,
        default=0.1,
        help="max star-match radius in arcsec (default 0.1 ~ 3 px)",
    )
    p_stack.add_argument("--threshold-sigma", type=float, default=5.0)
    p_stack.add_argument("--bg-sigma", type=float, default=3.0)
    p_stack.add_argument(
        "--interp-order",
        type=int,
        default=3,
        help="spline order for applying the shift (3 = cubic, keeps the PSF sharp)",
    )
    p_stack.add_argument("--no-preview", action="store_true")

    p_cmp = sub.add_parser("compare", help="Compare a stack to an official i2d mosaic.")
    p_cmp.add_argument("--stack", required=True)
    p_cmp.add_argument("--i2d", required=True)
    p_cmp.add_argument("--outdir", default="out")
    p_cmp.add_argument("--star-radius", type=float, default=8.0)
    p_cmp.add_argument(
        "--tiled",
        action="store_true",
        help="stream tile by tile (required for the full Stage 4 mosaic)",
    )
    p_cmp.add_argument("--tile-px", type=int, default=1024)
    p_cmp.add_argument("--halo-px", type=int, default=24)
    p_cmp.add_argument(
        "--write-diff",
        action="store_true",
        help="also write the full difference FITS (~1.4 GB)",
    )
    p_cmp.add_argument("--max-stars", type=int, default=6000)

    p_mos = sub.add_parser(
        "mosaic", help="Build the full mosaic on the fixed grid (Stage 4)."
    )
    _add_common(p_mos)
    _add_selection(p_mos)
    p_mos.add_argument("--grid", default="out/grid.fits")
    p_mos.add_argument("--out", default="out/f200w_all_detectors_mosaic.fits")
    p_mos.add_argument("--sigma", type=float, default=3.0)
    p_mos.add_argument("--iterations", type=int, default=3)
    p_mos.add_argument("--tile-px", type=int, default=1024)
    p_mos.add_argument("--halo-px", type=int, default=24)
    p_mos.add_argument("--cache-frames", type=int, default=32)
    p_mos.add_argument(
        "--registration",
        default="out/mosaic_registration.json",
        help="reuse a previous registration solution instead of re-solving",
    )
    p_mos.add_argument(
        "--recompute-registration",
        action="store_true",
        help="ignore an existing registration JSON and re-solve",
    )
    p_mos.add_argument("--fwhm-px", type=float, default=3.0)
    p_mos.add_argument("--match-radius-arcsec", type=float, default=0.1)
    p_mos.add_argument("--threshold-sigma", type=float, default=5.0)
    p_mos.add_argument("--bg-sigma", type=float, default=3.0)
    p_mos.add_argument("--interp-order", type=int, default=3)
    p_mos.add_argument("--scratch-dir", default=None)
    p_mos.add_argument(
        "--no-register",
        action="store_true",
        help="pure-WCS mosaic (skip registration entirely)",
    )
    p_mos.add_argument(
        "--no-cross-visit",
        action="store_true",
        help="skip the cross-visit alignment stage (per-group registration only)",
    )
    p_mos.add_argument(
        "--gauge-visit",
        type=int,
        default=None,
        help=(
            "visit whose header WCS defines the mosaic's absolute frame. The "
            "cross-visit correction is relative, so this origin is a choice. "
            "Default is the lowest-numbered visit, which is a heuristic: set it "
            "explicitly after measuring the raw frame WCSs against this "
            "filter's official i2d, because a different filter can have a "
            "different visit displaced"
        ),
    )
    p_mos.add_argument(
        "--plan-only",
        action="store_true",
        help="report the resource envelope and exit without building",
    )

    p_gauge = sub.add_parser(
        "gauge",
        help=(
            "Step 2: measure the raw per-visit frame WCS offsets against an "
            "official i2d and report the anchor visit."
        ),
    )
    _add_common(p_gauge)
    _add_selection(p_gauge)
    p_gauge.add_argument(
        "--i2d",
        required=True,
        help="this filter's official combined *_i2d.fits",
    )
    p_gauge.add_argument(
        "--exposure",
        type=int,
        default=1,
        help="exposure number supplying the reference stars, per group (default 1)",
    )
    p_gauge.add_argument(
        "--match-radius-arcsec",
        type=float,
        default=gauge_mod.DEFAULT_MATCH_ARCSEC,
        help=(
            "max star-match radius in arcsec, converted at the i2d's measured "
            f"scale (default {gauge_mod.DEFAULT_MATCH_ARCSEC})"
        ),
    )
    p_gauge.add_argument("--fwhm-px", type=float, default=gauge_mod.DEFAULT_FWHM_PX)
    p_gauge.add_argument("--threshold-sigma", type=float, default=5.0)
    p_gauge.add_argument(
        "--min-matched", type=int, default=gauge_mod.DEFAULT_MIN_MATCHED
    )
    p_gauge.add_argument(
        "--max-i2d-stars", type=int, default=gauge_mod.DEFAULT_MAX_I2D_STARS
    )
    p_gauge.add_argument(
        "--max-frame-stars", type=int, default=gauge_mod.DEFAULT_MAX_FRAME_STARS
    )
    p_gauge.add_argument(
        "--json-out",
        default=None,
        help="per-group JSON (default <outdir>/<filter>_step2_detectors.json)",
    )
    p_gauge.add_argument(
        "--log-out",
        default=None,
        help="text report (default <outdir>/<filter>_step2.log)",
    )
    p_gauge.add_argument(
        "--no-write", action="store_true", help="print only, write no files"
    )

    p_dl = sub.add_parser("download", help="List and download MAST cal/i2d products.")
    p_dl.add_argument("--input-dir", default=download.DEFAULT_DATA_ROOT)
    p_dl.add_argument("--i2d-dir", default=download.DEFAULT_I2D_ROOT)
    p_dl.add_argument("--outdir", default="out")
    p_dl.add_argument("--manifest", default=None, help="manifest CSV path")
    p_dl.add_argument("--filters", default="", nargs="+", help="e.g. F200W")
    p_dl.add_argument("--detectors", default="", nargs="+", help="e.g. nrca1 nrcblong")
    p_dl.add_argument("--products", default="cal,i2d")
    p_dl.add_argument("--proposal-id", type=int, default=download.PROPOSAL_ID)
    p_dl.add_argument(
        "--stage2",
        action="store_true",
        help=(
            "pre-selected Stage-2 set: all cal files, the combined i2d mosaic "
            "and the per-exposure i2d validation subset, for each --filters "
            "filter (default: F200W)"
        ),
    )
    p_dl.add_argument(
        "--stage2-kinds",
        default="cal,combined_i2d,subset_i2d",
        help=(
            "comma-separated subset of the Stage-2 categories to plan and "
            "download: cal, combined_i2d, subset_i2d. The per-exposure i2d "
            "subset_i2d is only needed for per-exposure validation; Stage-4 "
            "mosaic validation uses combined_i2d alone."
        ),
    )
    p_dl.add_argument("--plan-only", action="store_true", help="print the plan, do not download")
    p_dl.add_argument("--yes", action="store_true", help="skip the interactive prompt")

    p_verify = sub.add_parser(
        "verify",
        help="Check the cal/i2d roots against the MAST plan and rebuild manifests.",
    )
    p_verify.add_argument("--input-dir", default=download.DEFAULT_DATA_ROOT)
    p_verify.add_argument("--i2d-dir", default=download.DEFAULT_I2D_ROOT)
    p_verify.add_argument("--outdir", default="out")
    p_verify.add_argument("--proposal-id", type=int, default=download.PROPOSAL_ID)
    p_verify.add_argument(
        "--filter",
        dest="filter_name",
        nargs="+",
        default=None,
        help=(
            "filters to verify (default: every NIRCam filter in the MAST plan). "
            "Must cover every filter present under the data root, or those "
            "files are reported as 'extra' and the exit code stops meaning "
            "anything"
        ),
    )
    p_verify.add_argument(
        "--require",
        nargs="+",
        default=None,
        metavar="FILTER",
        help=(
            "only these filters gate the exit code (e.g. --require F090W). "
            "Absent products for other filters are still listed in the report "
            "but do not fail the run; the manifests still cover every filter"
        ),
    )
    p_verify.add_argument(
        "--require-kind",
        nargs="+",
        default=None,
        metavar="KIND",
        help=(
            "only these product kinds gate the exit code (CAL, I2D); combine "
            "with --require for a partially downloaded filter"
        ),
    )
    p_verify.add_argument(
        "--no-rebuild",
        action="store_true",
        help="report only, do not rewrite the manifests",
    )

    p_grid = sub.add_parser(
        "grid", help="Build (once) or describe the fixed common output grid."
    )
    p_grid.add_argument("--outdir", default="out")
    p_grid.add_argument("--scale", type=float, default=grid_mod.DEFAULT_PIXEL_SCALE_ARCSEC)
    p_grid.add_argument("--pad", type=int, default=grid_mod.DEFAULT_PAD_PX)
    p_grid.add_argument("--proposal-id", type=int, default=download.PROPOSAL_ID)
    p_grid.add_argument("--refresh-cache", action="store_true")
    p_grid.add_argument("--show", action="store_true")

    return parser


def _load_exposures(input_dir: str) -> list[io.CalExposure]:
    files = io.find_cal_files(input_dir)
    if not files:
        sys.exit(f"no *_cal.fits files found under {input_dir}")
    return [io.read_cal_exposure(f) for f in files]


def _select_exposures(
    exposures: list[io.CalExposure],
    detector: list[str] | None,
    filter_name: list[str] | None,
    pupil: list[str] | None = None,
) -> list[io.CalExposure]:
    """Keep only the requested detectors / filters / pupils.

    ``FILTER`` is ``F444W`` for both the CLEAR and the F470N exposure, so
    ``--filter F444W`` alone used to plan and build one 80-frame mosaic from two
    different bandpasses without a word of complaint.  The check runs on the
    *result*, not on whether ``--pupil`` was passed, so asking for both pupils
    explicitly is refused too rather than honoured.
    """
    selected = io.select_exposures(
        exposures, detector=detector, filter_name=filter_name, pupil=pupil
    )
    found = io.mixed_pupils(selected)
    if len(found) > 1:
        raise SystemExit(
            f"the selection spans {len(found)} pupils ({', '.join(found)}), so it "
            f"would mix bandpasses in one stack: {', '.join(found)} share "
            f"FILTER=F444W and are only distinguishable by PUPIL. Pass "
            f"--pupil with exactly one of them."
        )
    return selected


def _load_selected(args: argparse.Namespace) -> list[io.CalExposure]:
    exposures = _select_exposures(
        _load_exposures(args.input_dir),
        args.detector,
        args.filter_name,
        getattr(args, "pupil", None),
    )
    if not exposures:
        sys.exit(
            "no exposures match the requested detector/filter/pupil "
            f"({args.detector}, {args.filter_name}, "
            f"{getattr(args, 'pupil', None)}) under {args.input_dir}"
        )
    detectors = sorted({e.detector for e in exposures})
    filters = sorted({e.filter_name for e in exposures})
    print(
        f"loaded {len(exposures)} exposures "
        f"({len(detectors)} detector(s): {', '.join(detectors)}; "
        f"filter(s): {', '.join(filters)}; pupil: "
        f"{', '.join(io.mixed_pupils(exposures))})"
    )
    return exposures


def run_gauge(args: argparse.Namespace) -> None:
    exposures = _load_selected(args)
    label = gauge_mod.filter_label(exposures)
    print(
        f"gauge target: {label} "
        f"({len(exposures)} exposures, reference exposure {args.exposure})"
    )
    try:
        result = gauge_mod.measure_gauge(
            exposures,
            args.i2d,
            reference_exposure=args.exposure,
            match_arcsec=args.match_radius_arcsec,
            fwhm_px=args.fwhm_px,
            threshold_sigma=args.threshold_sigma,
            max_i2d_stars=args.max_i2d_stars,
            max_frame_stars=args.max_frame_stars,
            min_matched=args.min_matched,
        )
    except gauge_mod.GaugeError as exc:
        sys.exit(f"gauge failed: {exc}")

    print()
    print(gauge_mod.format_gauge_report(result, filter_label=label))

    if args.no_write:
        return
    outdir = Path(args.outdir)
    # F444W;F470N must not land in the same file as F444W:CLEAR, so the pupil
    # stays in the stem -- the two are different measurements.
    stem = label.lower().replace(";", "_") or "gauge"
    json_out = (
        Path(args.json_out) if args.json_out else outdir / f"{stem}_step2_detectors.json"
    )
    log_out = Path(args.log_out) if args.log_out else outdir / f"{stem}_step2.log"
    gauge_mod.write_gauge_outputs(
        result, json_path=json_out, log_path=log_out, filter_label=label
    )
    print()
    print(f"wrote {json_out}")
    print(f"wrote {log_out}")
    print(
        f"next: re-run the mosaic with --gauge-visit {result.anchor_visit} "
        f"(this is a measured anchor, not the lowest-visit heuristic)"
    )


def run_inspect(args: argparse.Namespace) -> None:
    exposures = _load_selected(args)
    table = grouping.frame_table(exposures)
    print(grouping.format_frame_table(table))


def run_group(args: argparse.Namespace) -> None:
    exposures = _load_selected(args)
    groups, matrix = grouping.connected_groups(
        exposures, threshold=args.threshold, sample=args.sample
    )
    print(grouping.format_groups(exposures, groups, matrix, args.threshold))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    report = outdir / "overlap_groups.txt"
    report.write_text(
        grouping.format_groups(exposures, groups, matrix, args.threshold) + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {report}")


def run_stack(args: argparse.Namespace) -> None:
    exposures = _load_selected(args)
    if args.visit == "all":
        visits = sorted({e.visit for e in exposures})
    else:
        visits = [int(args.visit)]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for visit in visits:
        group = [e for e in exposures if e.visit == visit]
        if not group:
            sys.exit(f"no exposures for visit {visit}")
        base = group[0].name.replace("_cal.fits", "")
        label = f"visit {visit} ({len(group)} frames)"
        print(f"\n= stacking {label} with {args.method}, "
              f"sigma={args.sigma}, iters={args.iterations}, "
              f"combine={args.combine}")
        started = time.perf_counter()
        stacked, coverage, out_wcs = stack_mod.stack_visit(
            group,
            sigma=args.sigma,
            iterations=args.iterations,
            combine=args.combine,
            method=args.method,
            pixel_scale_arcsec=args.pixel_scale,
            register=not args.no_register,
            fwhm_px=args.fwhm_px,
            match_radius_arcsec=args.match_radius_arcsec,
            threshold_sigma=args.threshold_sigma,
            bg_sigma=args.bg_sigma,
            interp_order=args.interp_order,
            outdir=outdir,
        )
        elapsed = time.perf_counter() - started
        print(f"  done in {elapsed:.1f}s  output {stacked.shape} px, "
              f"max depth {int(coverage.max())}")

        fits_path = outdir / f"{base}_stack.fits"
        _write_stack_fits(
            fits_path, stacked, coverage, out_wcs, group, args
        )
        print(f"  wrote {fits_path}")

        if not args.no_preview:
            png_path = plotting.save_preview_png(
                stacked,
                outdir / f"{base}_stack.png",
                title=f"NGC 3324 F200W nrca1 - {label} stack",
            )
            print(f"  wrote {png_path}")


def _write_stack_fits(
    path: Path,
    stacked: np.ndarray,
    coverage: np.ndarray,
    out_wcs: WCS,
    group: list[io.CalExposure],
    args: argparse.Namespace,
) -> None:
    primary = fits.PrimaryHDU(data=stacked.astype(np.float32))
    header = primary.header
    header.update(out_wcs.to_header())
    first = group[0]
    header["BUNIT"] = ("MJy/sr", "input SCI unit, preserved through stacking")
    header["FILTER"] = first.filter_name
    header["DETECTOR"] = first.detector
    header["NCOMBINE"] = (len(group), "exposures combined")
    header["SIGMA"] = (args.sigma, "sigma-clip threshold")
    header["ITERS"] = (args.iterations, "sigma-clip iterations")
    header["CLIPMETH"] = (args.method, "reprojection method")
    header["EXPTOT"] = (sum(e.effexptm for e in group), "total exposure time s")
    history = f"Stacked [{first.filter_name}] {len(group)} exposures with reproject"
    header["HISTORY"] = history
    for e in group:
        header["HISTORY"] = f"  {e.name}"
    coverage_hdu = fits.ImageHDU(data=coverage, name="COVERAGE")
    hdul = fits.HDUList([primary, coverage_hdu])
    hdul.writeto(path, overwrite=True)


def run_mosaic(args: argparse.Namespace) -> None:
    from jwst_stack import mosaic as mosaic_mod

    exposures = _load_selected(args)
    print(f"loaded {len(exposures)} exposures")

    grid_path = Path(args.grid)
    if not grid_path.exists():
        sys.exit(f"fixed grid not found: {grid_path} (run `grid` first)")

    plan = mosaic_mod.plan_mosaic(
        exposures,
        grid_path,
        tile_px=args.tile_px,
        halo_px=args.halo_px,
        cache_frames=args.cache_frames,
        scratch_dir=args.scratch_dir,
    )
    print(plan.text())
    needed = plan.scratch_gb + plan.output_gb
    if plan.free_disk_gb < needed * 1.15:
        sys.exit(
            f"insufficient disk: need ~{needed * 1.15:.1f} GB, "
            f"have {plan.free_disk_gb:.1f} GB"
        )
    if args.plan_only:
        print("\n--plan-only: nothing written")
        return

    registrations: dict = {}
    reg_path = Path(args.registration)
    reuse = reg_path.exists() and not args.recompute_registration
    if reuse and args.gauge_visit is not None:
        sys.exit(
            f"--gauge-visit {args.gauge_visit} cannot be applied to the cached "
            f"registration in {reg_path}, which already encodes its own gauge. "
            "Re-run with --recompute-registration to re-solve against the new "
            "anchor; silently reusing the old gauge would build a mosaic in a "
            "frame you did not ask for."
        )
    if args.no_register:
        print("\nregistration disabled (--no-register): pure WCS positions")
    elif reuse:
        registrations = mosaic_mod.read_registration(reg_path)
        print(
            f"\nreusing registration from {reg_path} "
            f"({len(registrations)} frames, "
            f"{sum(1 for r in registrations.values() if r.applied)} shifted)"
        )
    else:
        print("\nsolving registration on mosaic-grid crops (one group per visit x detector)")
        grid_wcs, grid_shape = mosaic_mod.load_grid(grid_path)
        solutions, summary = mosaic_mod.solve_registration(
            exposures,
            grid_wcs,
            grid_shape,
            fwhm_px=args.fwhm_px,
            match_radius_arcsec=args.match_radius_arcsec,
            threshold_sigma=args.threshold_sigma,
            bg_sigma=args.bg_sigma,
            interp_order=args.interp_order,
            cross_visit=not args.no_cross_visit,
            gauge_visit=args.gauge_visit,
        )
        mosaic_mod.write_registration(solutions, summary, reg_path)
        registrations = {s.filename: s for s in solutions}
        moved = [s for s in solutions if not s.is_reference and s.applied]
        rms = [s.residual_rms_px for s in moved if np.isfinite(s.residual_rms_px)]
        print(
            f"  wrote {reg_path}: {len(moved)}/{len(solutions) - 1} frames shifted, "
            f"median |shift| {np.median([np.hypot(s.dx, s.dy) for s in moved]) if moved else 0:.4f} px, "
            f"mean residual rms {np.mean(rms) if rms else float('nan'):.3f} px"
        )

    print()
    result = mosaic_mod.build_mosaic(
        exposures,
        grid_path,
        args.out,
        registrations=registrations,
        sigma=args.sigma,
        iterations=args.iterations,
        interp_order=args.interp_order,
        tile_px=args.tile_px,
        halo_px=args.halo_px,
        cache_frames=args.cache_frames,
        scratch_dir=args.scratch_dir,
    )
    print()
    print(result.text())


def run_compare(args: argparse.Namespace) -> None:
    if args.tiled:
        from jwst_stack import validation

        result = validation.validate_mosaic(
            args.stack,
            args.i2d,
            outdir=args.outdir,
            tile_px=args.tile_px,
            halo_px=args.halo_px,
            star_radius_px=args.star_radius,
            max_stars=args.max_stars,
            write_diff=args.write_diff,
        )
        print()
        print(result.text())
        return
    with fits.open(args.stack) as hdul:
        mine = hdul[0].data
        out_wcs = WCS(hdul[0].header)
    official = plotting.reproject_i2d(args.i2d, out_wcs)
    metrics = plotting.compare_stacks(mine, official, args.star_radius)
    print("agreement metrics (stack minus official, on overlapping pixels):")
    print(f"  overlap pixels      : {metrics.overlap_pixels}")
    print(f"  median difference   : {metrics.median_diff:.4g} MJy/sr")
    print(f"  MAD of difference   : {metrics.mad_diff:.4g} MJy/sr")
    print(f"  RMS of difference   : {metrics.rms_diff:.4g} MJy/sr")
    print(f"  matched stars       : {metrics.matched_stars}")
    print(f"  median star offset  : {metrics.median_star_offset_px:.2f} px")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    diff_img = mine - official
    diff_path = outdir / f"{Path(args.stack).stem}_vs_{Path(args.i2d).stem}"
    fits.HDUList(
        [
            fits.PrimaryHDU(data=diff_img.astype(np.float32), header=out_wcs.to_header()),
            fits.ImageHDU(data=official.astype(np.float32), name="OFFICIAL"),
        ]
    ).writeto(diff_path.with_suffix(".fits"), overwrite=True)
    plotting.save_difference_png(
        mine,
        official,
        diff_path.with_suffix(".png"),
        title="stack minus official i2d mosaic",
    )
    print(f"wrote {diff_path.with_suffix('.fits')} and .png")


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "inspect":
        run_inspect(args)
    elif args.command == "group":
        run_group(args)
    elif args.command == "stack":
        run_stack(args)
    elif args.command == "compare":
        run_compare(args)
    elif args.command == "mosaic":
        run_mosaic(args)
    elif args.command == "gauge":
        run_gauge(args)
    elif args.command == "download":
        if getattr(args, "stage2", False):
            # propagate the exit code: run_stage2 returns 1 on a failed disk
            # check, and silently exiting 0 there hides the failure
            raise SystemExit(download.run_stage2(args))
        else:
            raise SystemExit(download.run_download(args))
    elif args.command == "verify":
        raise SystemExit(download.run_verify(args))
    elif args.command == "grid":
        grid_mod.run_grid(args)


if __name__ == "__main__":
    main()