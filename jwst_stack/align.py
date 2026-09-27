"""Common output-grid construction and WCS-based reprojection."""

from __future__ import annotations

import numpy as np
from astropy.wcs import WCS
from reproject import reproject_exact, reproject_interp

from jwst_stack.grouping import footprint_polygon


def _plain_tan(ctype: str) -> str:
    """Drop any ``-SIP`` suffix from a CTYPE string.

    The output grid is a plain tangent-plane projection; it carries no SIP
    polynomial coefficients, so it must not advertise an SIP distortion.
    """
    if ctype.endswith("-SIP"):
        return ctype[:-4]
    return ctype


def build_output_wcs(
    exposures,
    margin_px: int = 32,
    pixel_scale_arcsec: float | None = None,
) -> tuple[WCS, tuple[int, int]]:
    """Build a plain-TAN grid covering the union of the frame footprints.

    The output grid shares the reference exposure's projection, CD matrix
    (hence orientation and native pixel scale) and CRVAL; only ``CRPIX`` is
    shifted so that the union of all footprints, plus ``margin_px`` pixels
    on every side, fits exactly.  Output pixel (i, j) corresponds to
    reference-frame pixel (i + dx, j + dy) with integer (dx, dy), so source
    frames that already share the reference orientation land on the grid
    without cumulative re-sampling error.

    If ``pixel_scale_arcsec`` is given the CD matrix is rescaled to that
    value (uniform in both axes) while preserving orientation.
    """
    reference = exposures[0]
    ref_wcs = reference.wcs

    min_x = min_y = np.inf
    max_x = max_y = -np.inf
    for exp in exposures:
        corners = footprint_polygon(exp.wcs, axes=exp.shape[::-1])
        px, py = ref_wcs.all_world2pix(corners[:, 0], corners[:, 1], 0)
        min_x = min(min_x, float(np.min(px)))
        max_x = max(max_x, float(np.max(px)))
        min_y = min(min_y, float(np.min(py)))
        max_y = max(max_y, float(np.max(py)))

    nx = int(np.ceil(max_x - min_x)) + 2 * margin_px
    ny = int(np.ceil(max_y - min_y)) + 2 * margin_px

    offset_x = min_x - margin_px
    offset_y = min_y - margin_px

    out_wcs = WCS(naxis=2)
    out_wcs.wcs.ctype = [_plain_tan(ref_wcs.wcs.ctype[0]), _plain_tan(ref_wcs.wcs.ctype[1])]
    out_wcs.wcs.cunit = [ref_wcs.wcs.cunit[0], ref_wcs.wcs.cunit[1]]
    out_wcs.wcs.crval = ref_wcs.wcs.crval
    out_wcs.wcs.crpix = ref_wcs.wcs.crpix - np.array([offset_x, offset_y])
    out_wcs.wcs.cd = ref_wcs.wcs.cd.copy()

    if pixel_scale_arcsec is not None:
        current = float(
            out_wcs.proj_plane_pixel_scales()[0].to_value("arcsec")
        )
        factor = current / pixel_scale_arcsec
        out_wcs.wcs.cd = out_wcs.wcs.cd * factor

    out_wcs.array_shape = (ny, nx)
    return out_wcs, (ny, nx)


def reproject_to_grid(
    data: np.ndarray, input_wcs: WCS, output_wcs: WCS, method: str = "interp"
) -> tuple[np.ndarray, np.ndarray]:
    """Reproject *data* onto the common output grid.

    Returns (reprojected, footprint) where footprint[i, j] == 1 where the
    output pixel is inside the source frame and 0 elsewhere.  Every output
    pixel outside the footprint is forced to NaN so uncovered pixels are
    never treated as valid science in the stack.
    """
    shape_out = tuple(int(s) for s in output_wcs.array_shape)
    if method == "exact":
        result, footprint = reproject_exact(
            (data, input_wcs),
            output_projection=output_wcs,
            shape_out=shape_out,
            return_footprint=True,
        )
    else:
        result, footprint = reproject_interp(
            (data, input_wcs),
            output_projection=output_wcs,
            shape_out=shape_out,
            return_footprint=True,
        )
    result = np.array(result, copy=True)
    result[footprint <= 0] = np.nan
    return result, footprint