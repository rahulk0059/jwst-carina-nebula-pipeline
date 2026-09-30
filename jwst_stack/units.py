"""Explicit angle and pixel-unit conversions, with the pixel basis named.

This module exists because "0.07 px" is not a number.  A pixel offset only
becomes an angle once you multiply it by a scale, and this project has **three
different pixel frames** in play at once, whose scales differ by a factor of
~2.03:

``i2d_px``
    Pixels of an official combined ``*_i2d.fits``.  ``compare``'s star-match
    offsets live here, because the match runs in i2d pixel space.  The scale is
    measured per filter and is *not* a round number: 0.031227 for the
    short-wave bandpasses, and 0.062904 (F335M) / 0.062908 (F444W;CLEAR) /
    0.062936 (F444W;F470N) for the long-wave ones.

``grid_px``
    Pixels of the output grid the mosaic is written on.  ``mosaic``'s
    registration shifts and residual rms live here, because the solve runs on
    grid crops.  The scale is 0.031 (``out/grid.fits``) or 0.0629
    (``out/grid_f335m/grid.fits``).

``degrees``
    FITS ``CDi_j``/``CRVAL`` values and MAST ``s_region`` polygons, which are
    in degrees.  Degrees to arcsec is x**3600**.

Using the wrong scale is a silent factor-of-2.029 error, and it has already
happened here: the long-wave 0.031-grid cross-checks were once reported as
"agreeing to 11%" when they had been converted with the *grid* scale instead of
the *i2d* scale, understating each offset by 50.7% and manufacturing agreement
between two numbers that did not agree.  The same session also multiplied an
angle in degrees by 1e6 instead of 3600 and turned a 0.24" disagreement into a
bogus 67.7", which looked alarming enough to threaten a decision that was
actually fine.

So the rule this module enforces is narrow and deliberate: :func:`px_to_arcsec`
**requires** a basis and refuses an unknown one, and :func:`format_px` refuses
to print a pixel quantity without one.  There is deliberately no default.  A
conversion with an unstated basis is the exact failure this module is for, and
making it easy to call by accident is how it happened.

See "Gotchas when re-measuring astrometry here" in AGENTS.md.
"""

from __future__ import annotations

from typing import Final

import numpy as np

#: Degrees to arcsec.  NOT 1e6 - microarcsec, and the other half of that bug.
ARCSEC_PER_DEG: Final[float] = 3600.0

#: Basis name for pixels of an official combined ``*_i2d.fits``.
I2D_PX: Final[str] = "i2d_px"

#: Basis name for pixels of the mosaic output grid.
GRID_PX: Final[str] = "grid_px"

#: The pixel frames this project actually has.  Closed on purpose: an unrecognised
#: basis is far more likely to be a typo than a genuinely new frame.
PIXEL_BASES: Final[tuple[str, ...]] = (I2D_PX, GRID_PX)


def _check_basis(basis: str) -> str:
    if basis not in PIXEL_BASES:
        raise ValueError(
            f"unknown pixel basis {basis!r}; expected one of "
            f"{', '.join(PIXEL_BASES)}. A pixel offset is not an angle until its "
            "basis is named - see jwst_stack.units and AGENTS.md."
        )
    return basis


def deg_to_arcsec(degrees: float | np.ndarray) -> float | np.ndarray:
    """Convert degrees to arcsec (x3600, not x1e6)."""
    return degrees * ARCSEC_PER_DEG


def arcsec_to_deg(arcsec: float | np.ndarray) -> float | np.ndarray:
    """Convert arcsec to degrees."""
    return arcsec / ARCSEC_PER_DEG


def px_to_arcsec(
    value: float | np.ndarray,
    scale_arcsec_per_px: float,
    basis: str,
) -> float | np.ndarray:
    """Convert a pixel offset to arcsec, given the scale of *its own* basis.

    *basis* has no default and is validated, because the only two plausible
    values here differ by 2.03x and picking the wrong one is silent.
    """
    _check_basis(basis)
    if not np.isfinite(scale_arcsec_per_px) or scale_arcsec_per_px <= 0:
        raise ValueError(
            f"pixel scale must be positive and finite, got {scale_arcsec_per_px!r}"
        )
    return value * scale_arcsec_per_px


def px_label(basis: str) -> str:
    """Short human label for a basis, e.g. ``i2d px``."""
    _check_basis(basis)
    return basis.replace("_", " ")


def format_px(
    value: float,
    scale_arcsec_per_px: float,
    basis: str,
    digits: int = 4,
    arcsec_digits: int = 5,
) -> str:
    """Render a pixel offset with its basis and its arcsec equivalent.

    The output names the basis and the scale it used, so a number can never be
    read off a report without the information needed to interpret it::

        0.0818 i2d_px @ 0.062936 arcsec/px = 0.00515 arcsec
    """
    _check_basis(basis)
    if value is None or not np.isfinite(value):
        return "n/a"
    arcsec = px_to_arcsec(float(value), scale_arcsec_per_px, basis)
    return (
        f"{value:.{digits}f} {basis} @ {scale_arcsec_per_px:.6f} arcsec/px "
        f"= {arcsec:.{arcsec_digits}f} arcsec"
    )


def format_pair_px(
    dx: float,
    dy: float,
    scale_arcsec_per_px: float,
    basis: str,
    digits: int = 4,
    arcsec_digits: int = 5,
) -> str:
    """Render a ``(dx, dy)`` pair with basis and arcsec, for report lines."""
    _check_basis(basis)
    if not (np.isfinite(dx) and np.isfinite(dy)):
        return "n/a"
    ax, ay = px_to_arcsec(float(dx), scale_arcsec_per_px, basis), px_to_arcsec(
        float(dy), scale_arcsec_per_px, basis
    )
    return (
        f"dx {dx:+.{digits}f} dy {dy:+.{digits}f} {basis} @ "
        f"{scale_arcsec_per_px:.6f} arcsec/px = dx {ax:+.{arcsec_digits}f} "
        f"dy {ay:+.{arcsec_digits}f} arcsec"
    )
