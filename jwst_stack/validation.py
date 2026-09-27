"""Tiled validation of a Stage 4 mosaic against the official combined i2d.

The mosaic and the official ``jw02731-o001_..._f200w_i2d.fits`` are both about
1.4 GB per float32 image on the 15895 x 22130 grid, so the comparison is also
streamed tile by tile.  For every tile the official image is reprojected onto
the mosaic grid, the difference is accumulated, and stars are detected
independently in both images.  The report separates three things:

* **background difference** - robust medians of each image over their overlap
  and the median of ``mine - official``;
* **residual noise** - MAD and RMS of the difference, both over the whole
  overlap and restricted to source-free pixels, so stars, cosmic rays and
  nebulosity do not inflate the quoted "noise";
* **star position offsets** - the signed median ``dx``/``dy`` of matched stars
  plus the residual scatter, which is the number that actually tests whether
  the mosaic grid and the registration put stars in the right place.

Per-tile robust statistics are used for the source-free mask rather than one
global threshold, because NGC 3324 sits on a bright nebular gradient and a
single global cut would classify half the field as "source".
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from scipy.spatial import cKDTree

from jwst_stack.mosaic import (
    _reproject,
    free_disk_gb,
    iter_tiles,
    peak_rss_gb,
    tile_sub_wcs,
)
from jwst_stack.register import detect_stars

DEFAULT_TILE_PX = 1024
DEFAULT_HALO_PX = 24
DEFAULT_STAR_RADIUS_PX = 8.0
DEFAULT_MAX_STARS = 6000
_CLEAN_SIGMA = 3.0
_VALUE_CAP = 200_000_000


@dataclass
class ValidationResult:
    """All Stage 4 comparison metrics."""

    mosaic_path: str
    i2d_path: str
    shape: tuple[int, int]
    i2d_shape: tuple[int, int]
    overlap_pixels: int
    mosaic_fraction: float
    i2d_fraction: float
    n_tiles: int
    background_mosaic: float = float("nan")
    background_i2d: float = float("nan")
    background_difference: float = float("nan")
    diff_median: float = float("nan")
    diff_mad: float = float("nan")
    diff_rms: float = float("nan")
    diff_mean: float = float("nan")
    diff_p01: float = float("nan")
    diff_p99: float = float("nan")
    clean_pixels: int = 0
    clean_fraction: float = 0.0
    clean_median: float = float("nan")
    clean_mad: float = float("nan")
    clean_rms: float = float("nan")
    stars_mosaic: int = 0
    stars_i2d: int = 0
    matched_stars: int = 0
    star_dx_median: float = float("nan")
    star_dy_median: float = float("nan")
    star_dx_rms: float = float("nan")
    star_dy_rms: float = float("nan")
    star_offset_median: float = float("nan")
    star_offset_p16: float = float("nan")
    star_offset_p84: float = float("nan")
    star_offset_mad: float = float("nan")
    star_nominal_frac: float = float("nan")
    seconds: float = 0.0
    peak_rss_gb: float = 0.0
    scratch_gb: float = 0.0
    exact_values: int = 0
    subsampled: bool = False
    diff_path: str | None = None
    png_path: str | None = None

    def as_dict(self) -> dict:
        out = dict(self.__dict__)
        out["shape"] = list(self.shape)
        out["i2d_shape"] = list(self.i2d_shape)
        return out

    def text(self) -> str:
        def f(v, unit="", digits=4):
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                return "n/a"
            return f"{v:.{digits}f}{unit}"

        lines = [
            "Stage 4 validation: mosaic vs official combined i2d",
            f"  mosaic            : {self.mosaic_path}",
            f"  official          : {self.i2d_path}",
            f"  grid shape        : {self.shape[0]} x {self.shape[1]}",
            f"  i2d native shape  : {self.i2d_shape[0]} x {self.i2d_shape[1]}",
            f"  tiles             : {self.n_tiles}",
            f"  overlap pixels    : {self.overlap_pixels:,} "
            f"({self.mosaic_fraction:.1f}% of mosaic, {self.i2d_fraction:.1f}% of i2d)",
            "",
            "background",
            f"  median mosaic     : {f(self.background_mosaic, ' MJy/sr')}",
            f"  median official   : {f(self.background_i2d, ' MJy/sr')}",
            f"  median difference : {f(self.background_difference, ' MJy/sr')}",
            "",
            "residual noise (whole overlap)",
            f"  MAD               : {f(self.diff_mad, ' MJy/sr')}",
            f"  RMS               : {f(self.diff_rms, ' MJy/sr')}",
            f"  mean              : {f(self.diff_mean, ' MJy/sr')}",
            f"  1st / 99th pct    : {f(self.diff_p01)} / {f(self.diff_p99)}",
            "",
            f"residual noise (source-free pixels, {_CLEAN_SIGMA:.0f}-sigma per tile)",
            f"  clean pixels      : {self.clean_pixels:,} ({100 * self.clean_fraction:.1f}%)",
            f"  median            : {f(self.clean_median, ' MJy/sr')}",
            f"  MAD               : {f(self.clean_mad, ' MJy/sr')}",
            f"  RMS               : {f(self.clean_rms, ' MJy/sr')}",
            "",
            "star position offsets",
            f"  stars mine / i2d  : {self.stars_mosaic} / {self.stars_i2d}",
            f"  matched           : {self.matched_stars}",
            f"  median dx, dy     : {f(self.star_dx_median, ' px')} , "
            f"{f(self.star_dy_median, ' px')}",
            f"  median |offset|   : {f(self.star_offset_median, ' px')}",
            f"  offset MAD        : {f(self.star_offset_mad, ' px')}",
            f"  offset 16-84 pct  : {f(self.star_offset_p16)} / {f(self.star_offset_p84)}",
            f"  rms dx, dy        : {f(self.star_dx_rms, ' px')} , "
            f"{f(self.star_dy_rms, ' px')}",
            f"  within 0.5 px     : {100 * self.star_nominal_frac:.1f}%",
            "",
            f"  diff values used  : {self.exact_values:,}"
            + (" (subsampled)" if self.subsampled else " (exact)"),
            f"  wall time         : {self.seconds / 60.0:.1f} min",
            f"  peak RSS          : {self.peak_rss_gb:.2f} GB",
        ]
        if self.scratch_gb:
            lines.append(f"  scratch (peak)    : {self.scratch_gb:.2f} GB")
        if self.diff_path:
            lines.append(f"  wrote difference  : {self.diff_path}")
        if self.png_path:
            lines.append(f"  wrote preview     : {self.png_path}")
        return "\n".join(lines)


def _sigma_mad(values: np.ndarray) -> float:
    if values.size == 0:
        return float("nan")
    med = float(np.median(values))
    return 1.4826 * float(np.median(np.abs(values - med)))


def _open_i2d(path: str | Path):
    """Return (SCI array, WCS, (ny, nx)) for an official i2d mosaic.

    The science image is taken from the first HDU carrying 2-D data, so both
    the usual ``SCI`` extension and a science image in the primary header work.
    """
    with fits.open(path, memmap=True) as hdul:
        hdu = None
        for candidate in hdul:
            data = candidate.data
            if data is not None and getattr(data, "ndim", 0) == 2:
                hdu = candidate
                break
        if hdu is None:
            raise ValueError(f"no 2-D image extension found in {path}")
        header = hdu.header
        data = hdu.data
        wcs = WCS(header)
        shape = (int(header["NAXIS2"]), int(header["NAXIS1"]))
    return data, wcs, shape


def _tile_stars(
    mine_exp: np.ndarray, off_exp: np.ndarray, y0, y1, x0, x1, ey0, ex0, max_stars
):
    """Detect stars in the halo-expanded tile, keeping only interior sources.

    Returns two ``(x, y, flux)`` tuples in *mosaic-grid* pixel coordinates.
    The halo lets DAOStarFinder see whole stars near a tile edge; detections
    whose centroid falls outside the interior are discarded so a star is never
    counted twice.
    """
    mfin = np.isfinite(mine_exp)
    ofin = np.isfinite(off_exp)
    both = mfin & ofin
    if not np.any(both):
        empty = (np.empty(0), np.empty(0), np.empty(0))
        return empty, empty
    out = []
    for arr, fin in ((mine_exp, mfin), (off_exp, ofin)):
        data = np.where(both, arr, np.nan)
        xs, ys, fs = detect_stars(data, max_stars=max_stars)
        keep = (
            (xs >= (x0 - ex0))
            & (xs < (x1 - ex0))
            & (ys >= (y0 - ey0))
            & (ys < (y1 - ey0))
        )
        out.append((xs[keep] + ex0, ys[keep] + ey0, fs[keep]))
    return out[0], out[1]


def _match_tile_stars(xo, yo, xm, ym, radius_px: float) -> dict:
    """Match official stars to mosaic stars and summarise the signed offsets.

    Each official star is matched to its nearest mosaic star within
    *radius_px*; a mosaic star can be claimed by at most one official star, and
    matches are accepted nearest-first so a crowded pair cannot produce two
    offsets from the same neighbour.
    """
    out = {
        "n_matched": 0,
        "dx_median": float("nan"),
        "dy_median": float("nan"),
        "dx_rms": float("nan"),
        "dy_rms": float("nan"),
        "offset_median": float("nan"),
        "offset_p16": float("nan"),
        "offset_p84": float("nan"),
        "offset_mad": float("nan"),
        "nominal_frac": float("nan"),
    }
    xo = np.asarray(xo, dtype=float)
    yo = np.asarray(yo, dtype=float)
    xm = np.asarray(xm, dtype=float)
    ym = np.asarray(ym, dtype=float)
    if xo.size == 0 or xm.size == 0:
        return out
    tree = cKDTree(np.column_stack([xm, ym]))
    dist, idx = tree.query(np.column_stack([xo, yo]), k=1)
    order = np.argsort(dist)
    used: set[int] = set()
    dxs, dys = [], []
    for j in order:
        if dist[j] > radius_px:
            break
        i = int(idx[j])
        if i in used:
            continue
        used.add(i)
        dxs.append(float(xm[i] - xo[j]))
        dys.append(float(ym[i] - yo[j]))
    if not dxs:
        return out
    dx = np.asarray(dxs)
    dy = np.asarray(dys)
    off = np.hypot(dx, dy)
    out.update(
        n_matched=len(dx),
        dx_median=float(np.median(dx)),
        dy_median=float(np.median(dy)),
        dx_rms=float(np.sqrt(np.mean(dx**2))),
        dy_rms=float(np.sqrt(np.mean(dy**2))),
        offset_median=float(np.median(off)),
        offset_p16=float(np.percentile(off, 16)),
        offset_p84=float(np.percentile(off, 84)),
        offset_mad=_sigma_mad(off),
        nominal_frac=float(np.mean(off <= 0.5)),
    )
    return out


def _save_preview(diff_mm, out_png: Path, target: int = 1800) -> Path:
    """Block-reduce the difference memmap and write an asinh PNG."""
    ny, nx = diff_mm.shape
    fy = max(1, ny // target)
    fx = max(1, nx // target)
    rows = max(1, ny // fy)
    cols = max(1, nx // fx)
    small = np.full((rows, cols), np.nan, dtype=np.float64)
    for r in range(rows):
        y0 = r * fy
        y1 = min(ny, y0 + fy)
        block = np.asarray(diff_mm[y0:y1, : cols * fx], dtype=np.float64)
        if block.size == 0:
            continue
        bh = block.shape[0] // fy
        if bh == 0:
            small[r, : block.shape[1] // fx] = np.nanmean(
                block.reshape(block.shape[0], 1, -1, fx), axis=(0, 3)
            )
            continue
        usable = (block.shape[0] // fy) * fy
        with np.errstate(invalid="ignore"):
            small[r, : block.shape[1] // fx] = np.nanmean(
                block[:usable].reshape(bh, fy, -1, fx), axis=(1, 3)
            )
    finite = small[np.isfinite(small)]
    vmax = float(np.percentile(np.abs(finite), 99.0)) if finite.size else 1.0
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    disp = np.arcsinh(small / vmax)
    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(disp, origin="lower", cmap="RdBu_r", interpolation="nearest",
                   vmin=-3, vmax=3)
    ax.set_title("mosaic minus official i2d (arcsinh, block-reduced)")
    ax.set_xlabel("mosaic pixel x")
    ax.set_ylabel("mosaic pixel y")
    fig.colorbar(im, ax=ax, label="arcsinh(mosaic - i2d) / p99")
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return out_png


def validate_mosaic(
    mosaic_path: str | Path,
    i2d_path: str | Path,
    outdir: str | Path | None = None,
    tile_px: int = DEFAULT_TILE_PX,
    halo_px: int = DEFAULT_HALO_PX,
    star_radius_px: float = DEFAULT_STAR_RADIUS_PX,
    max_stars: int = DEFAULT_MAX_STARS,
    write_diff: bool = False,
    scratch_dir: str | Path | None = None,
    progress_every: int = 20,
) -> ValidationResult:
    """Compare a Stage 4 mosaic with the official i2d, tile by tile."""
    started = time.time()
    mosaic_path = Path(mosaic_path)
    outdir = Path(outdir) if outdir else mosaic_path.parent
    outdir.mkdir(parents=True, exist_ok=True)
    stem = f"{mosaic_path.stem}_vs_{Path(i2d_path).stem}"

    with fits.open(mosaic_path, memmap=True) as hdul:
        header = hdul[0].header
        grid_wcs = WCS(header)
        shape = (int(header["NAXIS2"]), int(header["NAXIS1"]))
        mine = hdul[0].data
    official, i2d_wcs, i2d_shape = _open_i2d(i2d_path)

    scratch = Path(scratch_dir) if scratch_dir else outdir
    scratch.mkdir(parents=True, exist_ok=True)
    raw_path = scratch / "validate_diff.f32"
    diff_mm = np.memmap(raw_path, dtype=np.float32, mode="w+", shape=shape)
    scratch_gb = shape[0] * shape[1] * 4 / 1e9

    tiles = iter_tiles(shape, tile_px)
    print(
        f"validate: mosaic {shape[0]}x{shape[1]} vs i2d {i2d_shape[0]}x{i2d_shape[1]}, "
        f"{len(tiles)} tiles of {tile_px}px"
    )
    print(
        f"  scratch {scratch_gb:.2f} GB, free disk {free_disk_gb(scratch):.1f} GB"
    )

    overlap = 0
    sum_diff = sum_sq = sum_abs = 0.0
    sum_mine = sum_off = sum_mine_sq = sum_off_sq = 0.0
    n_finite = 0
    all_chunks: list[np.ndarray] = []
    clean_chunks: list[np.ndarray] = []
    n_all = 0
    n_clean = 0
    mx, my, mf = [], [], []
    ox, oy, of = [], [], []

    for ti, (y0, y1, x0, x1) in enumerate(tiles, start=1):
        ey0 = max(0, y0 - halo_px)
        ex0 = max(0, x0 - halo_px)
        ey1 = min(shape[0], y1 + halo_px)
        ex1 = min(shape[1], x1 + halo_px)
        exp_wcs = tile_sub_wcs(grid_wcs, ey0, ey1, ex0, ex1)
        off_exp, _ = _reproject(official, i2d_wcs, exp_wcs)
        del exp_wcs
        mine_exp = np.asarray(mine[ey0:ey1, ex0:ex1], dtype=np.float32)

        mine_t = np.ascontiguousarray(mine_exp[y0 - ey0 : y1 - ey0, x0 - ex0 : x1 - ex0])
        off_t = np.ascontiguousarray(off_exp[y0 - ey0 : y1 - ey0, x0 - ex0 : x1 - ex0])
        both = np.isfinite(mine_t) & np.isfinite(off_t)
        nb = int(both.sum())
        overlap += nb
        diff_mm[y0:y1, x0:x1] = mine_t - off_t
        if nb:
            mv = mine_t[both].astype(np.float64)
            ov = off_t[both].astype(np.float64)
            dv = mv - ov
            sum_diff += float(dv.sum())
            sum_sq += float(np.square(dv).sum())
            sum_abs += float(np.abs(dv).sum())
            sum_mine += float(mv.sum())
            sum_off += float(ov.sum())
            sum_mine_sq += float(np.square(mv).sum())
            sum_off_sq += float(np.square(ov).sum())
            n_finite += nb

            step = 1 if (n_all + nb) <= _VALUE_CAP else max(1, (n_all + nb) // _VALUE_CAP)
            chunk = dv[::step].astype(np.float32)
            all_chunks.append(chunk)
            n_all += chunk.size

            med_m = float(np.median(mv))
            mad_m = 1.4826 * float(np.median(np.abs(mv - med_m)))
            med_o = float(np.median(ov))
            mad_o = 1.4826 * float(np.median(np.abs(ov - med_o)))
            thr_m = max(_CLEAN_SIGMA * mad_m, 1e-12)
            thr_o = max(_CLEAN_SIGMA * mad_o, 1e-12)
            clean = (np.abs(mv - med_m) < thr_m) & (np.abs(ov - med_o) < thr_o)
            ncv = int(clean.sum())
            if ncv:
                cstep = 1 if (n_clean + ncv) <= _VALUE_CAP else max(1, (n_clean + ncv) // _VALUE_CAP)
                cchunk = dv[clean][::cstep].astype(np.float32)
                clean_chunks.append(cchunk)
                n_clean += cchunk.size
        (mxs, mys, mfs), (oxs, oys, ofs) = _tile_stars(
            mine_exp, off_exp, y0, y1, x0, x1, ey0, ex0, max_stars
        )
        if mxs.size:
            mx.append(mxs)
            my.append(mys)
            mf.append(mfs)
        if oxs.size:
            ox.append(oxs)
            oy.append(oys)
            of.append(ofs)
        del mine_t, off_t, both, mine_exp, off_exp

        if progress_every and (ti % progress_every == 0 or ti == len(tiles)):
            el = time.time() - started
            rate = ti / el if el else 0.0
            print(
                f"  tile {ti:>4}/{len(tiles)}  {el / 60.0:5.1f} min  "
                f"{rate:5.2f} tile/s  overlap {overlap:,}",
                flush=True,
            )

    values = np.concatenate(all_chunks) if all_chunks else np.empty(0, np.float32)
    clean_values = np.concatenate(clean_chunks) if clean_chunks else np.empty(0, np.float32)
    del all_chunks, clean_chunks

    star_x_m = np.concatenate(mx) if mx else np.empty(0)
    star_y_m = np.concatenate(my) if my else np.empty(0)
    star_f_m = np.concatenate(mf) if mf else np.empty(0)
    star_x_o = np.concatenate(ox) if ox else np.empty(0)
    star_y_o = np.concatenate(oy) if oy else np.empty(0)
    star_f_o = np.concatenate(of) if of else np.empty(0)
    if star_f_m.size > max_stars:
        keep = np.argsort(star_f_m)[::-1][:max_stars]
        star_x_m, star_y_m = star_x_m[keep], star_y_m[keep]
    if star_f_o.size > max_stars:
        keep = np.argsort(star_f_o)[::-1][:max_stars]
        star_x_o, star_y_o = star_x_o[keep], star_y_o[keep]

    star_fields = _match_tile_stars(
        star_x_o, star_y_o, star_x_m, star_y_m, star_radius_px
    )

    med = float(np.median(values)) if values.size else float("nan")
    result = ValidationResult(
        mosaic_path=str(mosaic_path),
        i2d_path=str(i2d_path),
        shape=shape,
        i2d_shape=i2d_shape,
        overlap_pixels=overlap,
        mosaic_fraction=100.0 * overlap / (shape[0] * shape[1]),
        i2d_fraction=100.0 * overlap / (i2d_shape[0] * i2d_shape[1]),
        n_tiles=len(tiles),
        background_mosaic=sum_mine / n_finite if n_finite else float("nan"),
        background_i2d=sum_off / n_finite if n_finite else float("nan"),
        background_difference=med,
        diff_median=med,
        diff_mad=_sigma_mad(values),
        diff_rms=float(np.sqrt(np.mean(np.square(values)))) if values.size else float("nan"),
        diff_mean=float(np.mean(values)) if values.size else float("nan"),
        diff_p01=float(np.percentile(values, 1)) if values.size else float("nan"),
        diff_p99=float(np.percentile(values, 99)) if values.size else float("nan"),
        clean_pixels=n_clean,
        clean_fraction=(n_clean / n_finite) if n_finite else 0.0,
        clean_median=float(np.median(clean_values)) if clean_values.size else float("nan"),
        clean_mad=_sigma_mad(clean_values),
        clean_rms=(
            float(np.sqrt(np.mean(np.square(clean_values))))
            if clean_values.size
            else float("nan")
        ),
        stars_mosaic=int(star_x_m.size),
        stars_i2d=int(star_x_o.size),
        matched_stars=int(star_fields["n_matched"]),
        star_dx_median=float(star_fields["dx_median"]),
        star_dy_median=float(star_fields["dy_median"]),
        star_dx_rms=float(star_fields["dx_rms"]),
        star_dy_rms=float(star_fields["dy_rms"]),
        star_offset_median=float(star_fields["offset_median"]),
        star_offset_p16=float(star_fields["offset_p16"]),
        star_offset_p84=float(star_fields["offset_p84"]),
        star_offset_mad=float(star_fields["offset_mad"]),
        star_nominal_frac=float(star_fields["nominal_frac"]),
        seconds=time.time() - started,
        peak_rss_gb=peak_rss_gb(),
        scratch_gb=scratch_gb,
        exact_values=int(values.size),
        subsampled=bool(values.size < overlap),
    )
    del values, clean_values

    diff_out = None
    png_out = None
    if write_diff:
        diff_out = outdir / f"{stem}_diff.fits"
        print("  writing difference FITS ...", flush=True)
        with fits.HDUList(
            [
                fits.PrimaryHDU(data=diff_mm, header=grid_wcs.to_header(relax=True)),
                fits.ImageHDU(data=np.zeros(shape, dtype=np.uint8), name="MASK"),
            ]
        ) as hdul:
            hdul.writeto(diff_out, overwrite=True)
        png_out = _save_preview(diff_mm, outdir / f"{stem}_diff.png")
    diff_mm.flush()
    del diff_mm
    try:
        raw_path.unlink()
    except OSError:
        pass
    result.diff_path = str(diff_out) if diff_out else None
    result.png_path = str(png_out) if png_out else None
    result.peak_rss_gb = max(result.peak_rss_gb, peak_rss_gb())
    result.seconds = time.time() - started

    (outdir / f"{stem}.json").write_text(
        json.dumps(result.as_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result
