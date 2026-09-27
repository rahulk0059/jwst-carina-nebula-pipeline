"""Per-visit alignment and sigma-clipped combination of exposures."""

from __future__ import annotations

import numpy as np
from astropy.stats import SigmaClip, sigma_clip
from pathlib import Path

from jwst_stack.align import build_output_wcs, reproject_to_grid
from jwst_stack.io import CalExposure, mask_scaled_sci
from jwst_stack.register import (
    RegistrationResult,
    format_registration_summary,
    register_overlays,
    write_registration_json,
)

_MIN_VALID_FOR_CLIP = 3


def stack_visit(
    exposures: list[CalExposure],
    sigma: float = 3.0,
    iterations: int = 3,
    combine: str = "median",
    method: str = "interp",
    pixel_scale_arcsec: float | None = None,
    margin_px: int = 32,
    register: bool = True,
    fwhm_px: float = 3.0,
    match_radius_arcsec: float = 0.1,
    threshold_sigma: float = 5.0,
    bg_sigma: float = 3.0,
    interp_order: int = 3,
    outdir: str | Path | None = None,
) -> tuple[np.ndarray, np.ndarray, object]:
    """Align and combine the group of exposures onto a shared grid.

    Returns (stack, coverage, output_wcs) where *coverage* counts how many
    exposures contribute a valid value to each output pixel.  ``combine`` is
    ``"median"`` or ``"mean"``; masking and sigma-clipping are always
    NaN-aware.

    With ``register=True`` (the default) photutils star registration is run on
    top of the WCS reprojections (see :func:`jwst_stack.register.
    register_overlays`): per-frame translation (dx, dy) and additive
    background offsets are solved in output-grid pixels and applied before
    combination.  Pass ``register=False`` to reproduce the pure-WCS stack
    exactly.  When *outdir* is given and registration runs,
    ``registration.json`` and ``registration_report.txt`` are written there
    with the per-frame solutions and residuals.
    """
    out_wcs, shape_out = build_output_wcs(
        exposures,
        margin_px=margin_px,
        pixel_scale_arcsec=pixel_scale_arcsec,
    )

    overlays = []
    footprints = []
    for exp in exposures:
        exp.load()
        cleaned = mask_scaled_sci(exp.sci, exp.dq)
        overlay, footprint = reproject_to_grid(
            cleaned, exp.wcs, out_wcs, method=method
        )
        overlays.append(overlay)
        footprints.append(footprint)
        exp.release()

    if register and len(overlays) >= 2:
        scale = pixel_scale_arcsec
        if scale is None:
            scale = float(out_wcs.proj_plane_pixel_scales()[0].to_value("arcsec"))
        overlays, reg_result = register_overlays(
            overlays,
            footprints,
            exposures,
            out_wcs,
            pixel_scale_arcsec=scale,
            fwhm_px=fwhm_px,
            match_radius_arcsec=match_radius_arcsec,
            threshold_sigma=threshold_sigma,
            bg_sigma=bg_sigma,
            interp_order=interp_order,
        )
        if outdir is not None:
            reg_dir = Path(outdir)
            reg_dir.mkdir(parents=True, exist_ok=True)
            write_registration_json(reg_result, reg_dir / "registration.json")
            report = format_registration_summary(reg_result)
            (reg_dir / "registration_report.txt").write_text(
                report + "\n", encoding="utf-8"
            )
            print(report)

    stack, coverage = sigma_clip_stack(overlays, sigma, iterations, combine)
    return stack, coverage, out_wcs


def sigma_clip_stack(
    overlays: list[np.ndarray],
    sigma: float = 3.0,
    iterations: int = 3,
    combine: str = "median",
) -> tuple[np.ndarray, np.ndarray]:
    """Combine already-reprojected frames with robust sigma clipping.

    NaN pixels are excluded; a pixel is only clipped when at least
    ``_MIN_VALID_FOR_CLIP`` frames contribute.  With fewer valid frames the
    plain median/mean of the valid values is used instead of clipping, which
    keeps a 2- or 3-frame overlap behaving sensibly.  Returns (stack,
    coverage) where coverage is the number of finite inputs per pixel.
    """
    cube = np.stack(overlays, axis=0)
    valid = np.isfinite(cube)
    coverage = valid.sum(axis=0).astype(np.int16)
    clipable = coverage >= _MIN_VALID_FOR_CLIP

    masked = sigma_clip(
        cube,
        sigma=sigma,
        maxiters=iterations,
        axis=0,
        masked=True,
        copy=True,
        cenfunc="median",
        stdfunc="mad_std",
    )
    clipped = masked.filled(np.nan)

    with np.errstate(invalid="ignore"):
        if combine == "mean":
            combined = np.nanmean(clipped, axis=0)
        else:
            combined = np.nanmedian(clipped, axis=0)

    keep = ~clipable
    if np.any(keep):
        if combine == "mean":
            fallback = np.nanmean(cube, axis=0)
        else:
            fallback = np.nanmedian(cube, axis=0)
        combined[keep] = fallback[keep]

    return combined, coverage