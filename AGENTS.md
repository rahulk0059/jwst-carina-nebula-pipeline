# AGENTS.md

Working notes for this repository. Read this before changing code or touching
the data roots.

## Project

`jwst_stack` stacks JWST NIRCam calibrated exposures (`*_cal.fits`, Level 2)
into mosaics for program **2731** (NGC 3324, Carina Nebula), NIRCam
short-wave. Per-visit stacks cover one detector (4 visits x 5 dithers = 20
exposures); the full **8-detector x 4-visit** mosaic (160 exposures) is Stage 4
and is the shipped product.

The heavy `jwst` package is **not** required. WCS comes from the `SCI` header
(`RA---TAN-SIP`, order 4), which `reproject` consumes directly.

## Current status

**All three 160-frame NIRCam filters are built and validated. Nothing is in
progress.** This is the cold-start summary; the sections after it hold the
evidence, and the per-filter sections are the authority for their own numbers.

| filter | mosaic (2.11 GB each) | gauge, measured | median offset vs official i2d | within 0.5 px | matched stars |
|--------|----------------------|-----------------|-------------------------------|--------------|---------------|
| F200W | `out/f200w_all_detectors_mosaic.fits` | visit 1, 0.0101 px, `explicit` | **0.0730 px** | 99.0% | 5882 / 6000 |
| F090W | `out/f090w_all_detectors_mosaic.fits` | visit 1, 0.017 px, `explicit` | **0.1384 px** | 78.5% | 5428 / 6000 |
| F187N | `out/f187n_all_detectors_mosaic.fits` | visit 1, 0.0060 px, `explicit` | **0.0647 px** | 99.3% | 5908 / 6000 |

Shared properties of all three: one grid (`out/grid.fits`, ny 15895 x
nx 22130), 352 tiles at 1024 px, 120.2-120.3 Mpx covered (34.2% of the grid),
max depth 17, ~30-36 min build each, peak RSS 4.47-4.58 GB. Each registration
JSON records `gauge_visit: 1` with `gauge_source: explicit`, because the gauge
was **measured** against that filter's own combined i2d before solving.

Three things a newcomer must not get wrong:

- **Do not re-gauge to the i2d.** The correction is relative, its absolute
  origin is a choice, and each filter's anchor was measured independently
  (0.0101 / 0.017 / 0.0060 px). Importing per-detector offsets from the same
  i2d would drive the reported residuals to ~0 *by construction* and destroy
  the only evidence that the pointing error exists. See "The gauge is the
  part that is easy to get wrong".
- **~0.05 px is a floor, not noise.** F200W and F187N both leave 0.0508 and
  0.0485 px unexplained after subtracting their two independently measured
  terms - agreeing to 5% across filters differing 5x in brightness. It is
  uncharacterised, and it is the largest open term. See "The three filters
  together, and a ~0.05 px common floor".
- **The per-visit translation model cannot reach zero.** F090W's residual is
  dominated by a 0.1118 px per-detector step that is *unidentifiable* from
  internal data (0 same-visit and 0 same-detector edges in the overlap graph).
  That is a real limitation of the model, not a solve failure.

To re-check the shipped state from scratch:

```powershell
.venv\Scripts\python -m pytest -q tests                      # 129 pass
.venv\Scripts\python -m jwst_stack.cli verify --outdir out   # disk truth
.venv\Scripts\python -m jwst_stack.cli compare --tiled `
  --stack out\f187n_all_detectors_mosaic.fits `
  --i2d "C:\data\jwst_i2d\mastDownload\JWST\jw02731-o001_t017_nircam_clear-f187n\jw02731-o001_t017_nircam_clear-f187n_i2d.fits" `
  --outdir out                                               # ~4 min, expect 0.0647 px
```

## What's next

Nothing below is needed to use what is shipped. Pick one.

1. **Single-visit filters: F335M, F444W, F444W;F470N** (40 cal files each, a
   ~944 MB combined i2d each, none downloaded yet). This is a *different kind*
   of test, not more of the same:
   - **No gauge problem.** One visit means no inter-visit pointing error, so
     the cross-visit stage has nothing to solve and there is no anchor to
     choose. Expect the intra-detector field term to be the whole story, and
     do not expect the ~0.05 px floor to be the only thing left.
   - **New grid, and check the pixel scale.** These are not the same detector
     geometry as the F090W/F187N/F200W set, so `out/grid.fits` may not contain
     them. Do the step-1 `_grid_bbox` check first, and if it does not fit,
     **build a new grid at a different path** - rebuilding `out/grid.fits` in
     place would invalidate all three shipped mosaics. Keep both grids.
   - The long-wave detectors (`nrcalong`/`nrcblong`, currently absent from
     disk) are a further step again: different detector family, different
     pixel scale, so expect another grid decision.
2. **Rotation-aware per-group registration.** Scoped, deliberately **not**
   implemented - see "Deferred: rotation-aware per-group registration". Worth
   ~20% of F090W's residual and nothing of the unidentifiable inter-detector
   half. Do not start it mid-filter; the risks are already written down.
3. **Characterise the ~0.05 px floor.** It is now the largest unexplained term
   in the project and it is not physical noise. Candidates are the per-star
   centroid floor and PSF-difference bias against the drizzle; neither is
   characterised. This is measurement work, not a code fix.
4. **Pause for productization.** The code is at a natural stopping point: 129
   offline tests, no linter, three validated end-to-end results. If the goal
   becomes a reusable tool rather than a set of measurement results, the
   remaining work is packaging, a config file instead of eight `--detector`
   flags, and CI - not more astronomy.

## Environment

- venv: `.venv\Scripts\python` (Windows / PowerShell 7)
- Run everything through the venv interpreter, e.g.
  `.venv\Scripts\python -m jwst_stack.cli <command>`
- Tests: `.venv\Scripts\python -m pytest -q tests` (129 tests, all offline)
- There is **no** linter or type-checker configured (no ruff/mypy/flake8 in the
  venv, no `pyproject.toml`/`setup.cfg`). `ast.parse` + the test suite are the
  correctness gate.

## Commands

| command | purpose |
|---------|---------|
| `inspect` | per-file metadata table |
| `group` | footprint-overlap groups |
| `stack` | register + align + sigma-clip stack one visit or all |
| `compare` | compare a stack to an official `*_i2d.fits` |
| `download` | MAST product discovery and resumable download |
| `verify` | compare the data roots against the MAST plan, rebuild manifests |
| `grid` | build/describe the fixed common output grid |
| `mosaic` | Stage 4: all 8 detectors x 4 visits on the fixed grid, tiled |

`inspect`, `group` and `stack` accept `--detector` / `--filter`. **Always pass
them.** The cal root holds all 8 NIRCam detectors, so `--visit 1` without
`--detector nrca1 --filter F200W` silently groups 40 exposures from different
detectors onto one 4760x11526 grid (and takes ~15 min instead of ~60 s).

## Data roots and disk truth

- cal: `C:\data\jwst_cal\mastDownload\JWST`
- i2d: `C:\data\jwst_i2d\mastDownload\JWST`

Verified state as of the last `verify` run (F187N scoped gate, exit 0):

| root | files | size | status |
|------|-------|------|--------|
| cal | 480 | 56.435 GB (52.56 GiB) | 480 ok, 0 size_mismatch, 120 missing, 0 extra |
| i2d | 40 | 20.645 GB (19.23 GiB) | 40 ok, 0 size_mismatch, 566 missing, 0 extra |
| **total** | **520** | **77.080 GB** (71.79 GiB) | |

The 480 cal files are the short-wave set: 8 detectors x 60 exposures
(3 filters x 4 visits x 5 dithers) for F090W, F187N and F200W. The **120
"missing" cal files are the long-wave detectors** (`nrcalong`, `nrcblong`, 60
each) - deliberately never downloaded, since the shipped mosaics are
short-wave. Of the 40 i2d on disk, **3 are the combined all-visit products**
(the ones validation actually uses) and 37 are per-exposure F200W products
from the original validation subset. The 566 "missing" i2d are the per-exposure
products for the other filters plus the F335M/F444W combined products; they are
listed in the manifest on purpose, because the manifest is an inventory of what
MAST offers versus what is on disk, not a record of intent. Do not "fix" them
by downloading unless asked.

MAST's advertised sizes are not uniform, and `verify` requires an **exact**
per-product match: 440 of the 480 cal files are 117,573,120 B, while 40 F200W
files are 117,570,240 B. Both count as `ok`.

## Manifest rule

`out/stage2_cal_manifest.csv` and `out/stage2_i2d_manifest.csv` are **disk
truth, not intent**. They are generated by `verify` and must never be
hand-edited.

- After any download, deletion or move under either data root, run
  `.venv\Scripts\python -m jwst_stack.cli verify --outdir out`
  to re-scan and rebuild both manifests.
- The manifest `size` column holds the **actual on-disk byte count**, and
  `status` is one of `ok`, `size_mismatch`, `missing`, `extra`.
- A stale manifest must never be trusted over a fresh scan: earlier manifests
  listed 45 cal rows while 160 files existed on disk, and omitted 8 i2d files
  entirely, so a truncated file was invisible until `verify` existed.
- `verify` exits non-zero while any gated problem remains, so it can gate a
  pipeline. Use `--require FILTER` (and `--require-kind CAL|I2D`) to scope the
  **exit code** to the products you actually intended to fetch:

  ```
  .venv\Scripts\python -m jwst_stack.cli verify --require F090W --require-kind CAL --outdir out
  ```

  This is the normal invocation, because the MAST plan lists far more than any
  one filter's worth of files. Scoping the exit code does **not** narrow the
  scan: the manifests are always rebuilt over every filter, since they are
  whole-program disk truth. Absent out-of-scope products are still listed, under
  their own heading, so scoping never hides anything. `size_mismatch` and
  `extra` gate at any scope - a truncated or unexpected file is a real anomaly.
- `--require-kind` exists because a filter is often only *partially* downloaded:
  F090W has all 160 cal files but only the combined i2d, so `--require F090W`
  alone can never pass. Scoping the kind (`--require F090W --require-kind CAL`)
  lets the cal side gate while the deliberately-skipped per-exposure i2d
  products are listed and ignored.
- **Do not repeat `--require` to gate two filters.** It is `nargs="+"`, and
  argparse lets a later occurrence *overwrite* the earlier one, so
  `--require F187N --require F090W` silently gates F090W only. Pass one flag
  with both values: `--require F187N F090W`. Verified: the repeated form
  returns `['F090W']`, the single form returns `['F187N', 'F090W']`.

## Deferred: rotation-aware per-group registration

**Status: scoped, deliberately NOT implemented.** Logged here so it is not
rediscovered from scratch, and so nobody "fixes" the 0.1384 px by reaching for
it mid-way through the remaining filters.

**What:** extend the per-group (dither-to-dither) registration in
`solve_registration` from a fitted translation to translation + rotation/skew,
so the intra-detector field term can be absorbed instead of being averaged over.

**Why it is legitimate (unlike the per-detector term):** a field gradient is a
smooth function of sky position and is identifiable from internal star matches
alone - the 5 dithers of a group overlap the same sky, so each frame's own
catalogue supplies the positional lever arm. No external reference is needed, so
the validation against the i2d stays independent. The per-detector *step* cannot
be treated this way (rank 17 of 32; 0 same-visit and 0 same-detector edges).

**Expected gain:** F090W 0.1384 -> ~0.11 px, about 20%. The 0.1118 px
inter-detector term is untouched and is the floor, so this cannot reach zero.
F200W's total residual is already 0.0730 px, so there is little to win there.

**Not urgent because** the remaining gain is small, and the work is riskier than
the number suggests:

1. The term is currently measured **per frame**; correcting it needs the fit done
   per frame from that frame's own dither matches, not one number per group.
2. It overturns the standing "translation only, no rotation, skew or scale"
   decision under "Stacking and registration". That decision's cubic-shift mask
   handling was validated for translation specifically, and the NaN-poisoning
   gotcha above is the concrete risk: rotating a footprint with a masked edge
   re-exposes the exact failure that `test_nan_footprint_does_not_poison_cubic_shift`
   guards.
3. The sigma-clipping requirement is not optional. Per-star RMS is ~0.7 px against
   a MAD of 0.11 px, so an unclipped rotation fit is set by a handful of blended
   and saturated cores - which is what made the first attempt at measuring this
   term report a meaningless 0 significance.

If it is ever attempted, the diagnostic to re-run afterwards is
`within_detector`-style per-group trend fits, checking the field term drops from
0.0812 px toward the per-star MAD floor rather than merely changing sign.

## Gotchas that caused real bugs

- **Server may ignore `Range`.** If MAST answers `200` to a resume request the
  whole body is sent, so the partial must be truncated. Reopening with
  `seek(0)` on a handle opened `"ab"` does nothing (O_APPEND) and produces an
  oversized corrupt file. `download_url` now picks the mode *after* reading the
  response status. Regression test:
  `test_download_url_truncates_when_server_ignores_range`.
- **Cubic resampling propagates NaN.** `scipy.ndimage.shift` with `order >= 2`
  runs a global spline prefilter, so a single NaN outside a footprint makes the
  whole frame NaN. Registration fills the mask before shifting and restores it
  after. Regression test: `test_nan_footprint_does_not_poison_cubic_shift`.
- **Linear resampling blurs.** Applying a sub-pixel shift with `order=1` costs
  ~5% wider stellar FWHM on visit 1, cancelling the alignment gain. The default
  is `interp_order=3`.
- **DAOStarFinder columns** are `xcentroid` / `ycentroid`, not `x` / `y`.
- **Match radius is in arcsec but applied in pixels.** 1.0 arcsec at
  0.0311 arcsec/px is 32 px, which cross-matches unrelated stars in a crowded
  field (residual rms 3.4-5.5 px). The default is 0.1 arcsec (~3 px), giving
  residual rms 0.08-0.22 px.
- **Translation sign.** `_match_stars` returns `frame - reference`, so the
  correction applied to the frame is `(-dy, -dx)`.
- **`build_items` drops combined products.** A product whose name has no
  detector/exposure (the all-visit i2d mosaic) is skipped, because
  `parse_product_name` returns `None`. Use `build_i2d_items` for i2d, which
  includes the mosaic.
- **The combined mosaic is `jw02731-o001_...`, not `o017`.** Older notes
  recorded `o017`; MAST lists `jw02731-o001_t017_nircam_clear-f200w_i2d.fits`
  at 5,420,381,760 B.
- **`download_items` reported success on a failed transfer.** Every failure was
  printed per file but never counted, so `run_stage2`/`run_download` returned 0
  and printed "Stage-2 download complete" after 89 of 110 files had failed with
  DNS errors. A partial download looked finished. `item_failed()` now defines
  success as `ok`/`skipped` only, both entry points return 1 and print
  `INCOMPLETE ... Re-run the same command to resume`, and `cli.py` propagates
  the code. Regression tests:
  `test_item_failed_treats_only_ok_and_skipped_as_success`,
  `test_run_stage2_returns_nonzero_when_a_file_fails`,
  `test_run_download_returns_nonzero_when_a_file_fails`.
- **MAST drops DNS after a burst of ~20-45 files.** Fetching F187N, the first
  46 files arrived, then every subsequent `urlopen` failed with
  `[Errno 11001] getaddrinfo failed` until the process was re-run some minutes
  later; `mast.stsci.edu` resolved fine the whole time. It is transient, not a
  permissions or URL problem, and a single long run cannot finish. Retry the
  same command repeatedly - it is resumable and skips files already at the
  exact expected size, so passes accumulate. An interrupted transfer leaves a
  short partial on disk (e.g. 43,384,832 B of 117,573,120 B); the next pass
  truncates and re-fetches it rather than appending. Check
  `run_stage2`'s exit code, not its output wording, to decide whether to retry.
- **`--stage2` was hard-coded to F200W and ignored `--filters`.** Asking for
  `download --stage2 --filters F187N --plan-only` returned a plan labelled
  "F200W cal files", sized for F200W. The failure was silent and total: the
  plan matched **0** F187N products, so it would have reported
  `0 to download` and exited 0 having fetched nothing. `plan_stage2` now takes
  `filter_names` (one or several), `run_stage2` resolves it from `--filters`
  via `resolve_stage2_filters`, and the plan carries `plan["filters"]` so
  `format_stage2_plan` labels with the filter actually planned instead of a
  baked-in string. The old label also hard-coded the *count* ("all 160"), which
  is wrong for the 40-file single-visit filters. Regression tests:
  `test_plan_stage2_respects_filters`,
  `test_resolve_stage2_filters_prefers_requested`,
  `test_run_stage2_plan_output_labels_requested_filter`.
- **There are two `download` parsers, and a flag added to one is invisible to
  the other.** `download.parse_args` and `cli.build_parser()` each define their
  own `download` subparser. Adding `--stage2-kinds` to only `download.py` left
  the real CLI crashing with `ValueError: --stage2-kinds must name at least one
  category`, because `cli.py`'s Namespace had no such attribute. Any new
  download/verify flag must be added to **both** and asserted on both:
  `test_cli_download_parses_stage2_flags` and
  `test_run_stage2_tolerates_missing_stage2_kinds_attr` exist for this.
  This is the same trap as `--require`/`--require-kind` earlier.
- **`cli.py` discarded `run_stage2`'s return code.** A failed disk check
  returned 1 from `run_stage2` and the CLI still exited 0, so a download that
  refused to run looked successful. The `download` branch now raises
  `SystemExit` with the code, like `verify`.
- **`--stage2` plans a third category that is easy to forget.**
  `subset_i2d` is the 55-file per-exposure i2d validation subset (nrca1 for all
  four visits plus all eight detectors for visit 1), 6.53 GB for F187N. It is
  only needed for per-exposure validation; Stage-4 mosaic validation uses the
  combined i2d alone. `--stage2-kinds cal,combined_i2d` selects the 161-file
  scope actually used for F090W and F187N. The default is unchanged
  (`cal,combined_i2d,subset_i2d`), so existing invocations still plan all three.
- **The combined i2d needs `build_i2d_items`, not `--products i2d`.**
  `download --products i2d --plan-only` lists only the 160 per-exposure
  products and never the all-visit mosaic, because the plain selection path
  drops products whose name has no detector/exposure components. Use
  `--stage2` (now filter-aware) or `build_i2d_items` when you need the
  combined product.
- **Batched `read` calls can scramble file labels.** Read one file per call, or
  use `inspect` / `ast` from Python when the exact on-disk text matters.

## Stacking and registration

- Output grid: plain TAN, source pixel scale and orientation, 32-px margin over
  the union of the visit's footprints.
- Registration is **on by default** for `stack`; `--no-register` reproduces the
  pure-WCS result. Translation only - no rotation, skew or scale. **This is
  now known to be a real limitation rather than a safe default**: the F090W
  intra-detector field term below is exactly what a translation-only per-group
  fit cannot represent. A rotation-aware fit is scoped and deferred, not
  implemented - see "Deferred: rotation-aware per-group registration".
- Reference frame is the first / lowest exposure of the visit. A frame needs
  >= 3 matched stars to be shifted, otherwise it keeps its WCS position.
- Background handling is an additive, sigma-clipped median offset only - no
  slope or scale - to preserve the NGC 3324 nebular gradient.
- Each visit writes `registration.json` and `registration_report.txt` (dx, dy,
  n_detected, n_matched, residual rms/max, background offset) next to the
  stack.
- **Within a visit**, the measured median |shift| is 0.07-0.38 px, so the header
  WCS is already good to ~0.1 px and within-visit registration is FWHM-neutral
  (within 0.05%). Its value is making the residual pointing error measured and
  auditable. Do not claim it sharpens the stack.
- **That ~0.1 px figure is a WITHIN-VISIT result and must not be quoted as a
  dataset-wide astrometric accuracy.** It was originally written here as a
  general claim and was wrong: it came from per-visit shifts only, and nothing
  compared one visit to another. Measured across visits, the header WCS is off
  by **0.8115 px** (visit 1 versus the mean of visits 2/3/4; see
  "Inter-visit alignment" below). Per-visit stacks are unaffected because each
  `stack` run is single-visit; only the Stage 4 mosaic spans visits.

## Stage 4: tiled full mosaic

`mosaic` is **not** `stack` on a bigger canvas. A 160-overlay in-memory stack
would need ~225 GB, so `mosaic.py` streams 1024-px tiles through memmapped
scratch (SCI float32 + COVERAGE uint16, 2.11 GB) and writes one FITS at the end.

- Fixed grid `out/grid.fits`: **ny 15895 x nx 22130 = 351,756,350 px**, 0.0310
  arcsec/px plain TAN, 32-px margin over the union of all 160 footprints. It is
  **wider than it is tall**, and `load_grid` returns the shape as
  `(ny, nx)` - see the axis-order gotcha below before writing any code that
  indexes it.
- **Bilinear reprojection, not cubic.** `reproject`'s `spline_filter` runs
  globally for `order>=2` and a single DQ-masked NaN outside a footprint poisons
  the whole tile. `_reproject` uses `reproject_interp` with `order=1`, and only
  the *registration* shift uses cubic `ndi_shift` (mask filled before, restored
  after) - the same recipe validated for per-visit stacks.
- Registration is solved per `(visit, detector)` group (32 groups of 5) on
  small crops **of the mosaic grid**, so the stored shifts are already in
  mosaic-grid pixels and need no rescaling.
- **A second, cross-visit stage is required** (`solve_cross_visit_alignment`).
  The per-group solve above cannot see inter-visit pointing error at all, and
  on this dataset that error is 0.85 px. See "Inter-visit alignment" below.
  Pass `--no-cross-visit` to disable it and reproduce the buggy behaviour.
- `_grid_bbox` samples the 17x17 edge grid, but **TAN `all_world2pix` returns
  NaN** for a frame near the antipode; the function must filter non-finite
  corners (and reject when <3 remain) or `int(np.floor(nan))` raises.

### Measured Stage 4 run (160 frames, F200W, all 8 detectors)

```
.venv\Scripts\python -m jwst_stack.cli mosaic `
  --detector nrca1 nrca2 nrca3 nrca4 nrcb1 nrcb2 nrcb3 nrcb4 --filter F200W `
  --grid out\grid.fits --out out\f200w_all_detectors_mosaic.fits `
  --registration out\mosaic_registration.json --scratch-dir out `
  --recompute-registration
```

| quantity | value |
|----------|-------|
| frames combined | 160/160 |
| tiles | 352 at 1024 px (max 27 frames/tile, mean 5.5) |
| reproject calls | 1931 |
| covered pixels | 120,196,133 (34.2% of grid) |
| max depth | 17 |
| per-group registration | 128/159 frames shifted, median \|shift\| 0.1284 px, mean residual rms 0.156 px |
| cross-visit stage | 26 edges, residual rms 0.074 px, corrections v1 (0,0) v2 (+0.501,+0.544) v3 (+0.456,+0.737) v4 (+0.564,+0.750) px |
| wall time | ~36 min registration + 30.1 min build |
| peak RSS | 4.58 GB (planner predicted 2.47 GB - it under-counts) |
| scratch / output | 2.11 GB / 2.11 GB |
| frame loads | 491 (cache hits 1440) |

The planner's `estimated peak RSS` is a floor, not a forecast: the real peak is
~1.8x it. Quote measured peaks, never the planner figure.

### Stage 4 validation vs the official combined i2d

The official mosaic is **nested in a subdirectory**:
`C:\data\jwst_i2d\mastDownload\JWST\jw02731-o001_t017_nircam_clear-f200w\jw02731-o001_t017_nircam_clear-f200w_i2d.fits`
(5,420,381,760 B, native 8588 x 14342). `compare --tiled` reprojects it tile by
tile, so the 34.2%-filled mosaic is never loaded whole.

**These are the corrected numbers, after the cross-visit fix.** The
contaminated pre-fix values are kept in the comparison table below so the
regression stays visible.

| quantity | value |
|----------|-------|
| overlap | 120,190,998 px (34.2% of mosaic, **97.6% of i2d**) |
| background (median) | mosaic 2.7875 vs official 2.7506 -> **+0.0431 MJy/sr** |
| residual noise, whole overlap | MAD 0.0477, RMS 6.3639 MJy/sr |
| residual noise, source-free (92.0% of px) | median 0.0429, MAD 0.0452, **RMS 0.0873 MJy/sr** |
| stars mine/i2d | 6000 / 6000 (both hit the `max_stars` cap) |
| matched stars | 5882 |
| star offset (mine - official) | median dx **-0.0240**, dy **+0.0327** px; median \|offset\| **0.0730 px** |
| offset scatter | MAD 0.0430 px, rms dx 0.2438 / dy 0.2923 px |
| within 0.5 px | **99.0%** |

| metric | pre-fix | fix, consensus gauge | fix, visit-1 gauge |
|--------|---------|----------------------|-------------------|
| median dx, dy (px) | +0.4219, +0.6458 | +0.4832, +0.7124 | **-0.0240, +0.0327** |
| median \|offset\| (px) | 0.7938 | 0.8613 | **0.0730** |
| offset MAD (px) | 0.2275 | 0.0613 | **0.0430** |
| within 0.5 px | 27.9% | 0.2% | **99.0%** |
| whole-overlap RMS | 9.4365 | 11.9333 | **6.3639** |
| source-free RMS | 0.0954 | 0.0992 | **0.0873** |
| background difference | +0.0458 | +0.0458 | **+0.0431** |

Reading these honestly:

- The 0.79 px is **not** a grid/drizzle convention difference. It was our own
  inter-visit misalignment; see "Inter-visit alignment" below. Do not quote the
  pre-fix 0.79 px as an error budget.
- The two acceptance tests both pass: internal cross-visit residuals are 0.071
  (x) / 0.076 (y) px rms over the 26 edges, and the comparison against the
  official i2d is 0.0730 px median. The i2d therefore has no unresolved
  contribution at the 0.1 px level - the residual is measurement scatter.
- Whole-overlap RMS (6.36) is still dominated by outliers - saturated star
  cores and footprint edges - which is why rms dx/dy (0.24/0.29) is much larger
  than the 0.043 px MAD. The source-free RMS (0.087) is the meaningful noise
  figure and is consistent with the background MAD.
- The background difference (+0.043 MJy/sr) is **not** explained by the
  alignment fix; a rigid star translation cannot move a median background. It
  is a genuine photometric/calibration difference and is unchanged.
- The star sample is **capped at 6000 per side** (`--max-stars`), so it measures
  the brightest stars only; raise the cap for a population study.
- Validation is exact-valued, not subsampled: all 120,190,998 overlap pixels.

## Inter-visit alignment (the 0.81 px bug)

The 0.79 px "global offset" was accepted as a grid convention difference for one
round because it was plausible. It was wrong. A check showed it was our own
error.

**Root cause.** `solve_registration` registers each `(visit, detector)` group
against *its own* first exposure. Nothing ever compared one visit to another,
so inter-visit header-WCS pointing error passed straight into the mosaic.

**How it was measured.** Each group's reference frame was star-detected (677-
1562 stars each), converted to mosaic-grid pixels, and cross-matched pairwise.
Only **26 of 496** group pairs have >=20 matched stars, and **all 26 are
cross-visit - zero same-visit edges**, because within a visit the 8 detectors
tile disjoint sky. Widening the match radius 1 -> 10 px reproduces exactly the
same 26 edges, so this is not a selection-bias artefact.

| model | variance explained (x, y) | residual rms |
|--------|--------------------------|-------------|
| per-visit translation | **90.7%, 94.2%** | 0.079, 0.068 px |
| per-detector translation | -6.4%, 31.5% | 0.266, 0.234 px |

Fitted per-visit translations **relative to the consensus of v2/v3/v4**: v1
(-0.368, -0.485), v2 (+0.118, +0.037), v3 (+0.071, +0.228), v4 (+0.180, +0.220).
**Visit 1 is displaced 0.8115 px (dx -0.491, dy -0.646) from the mean of visits
2/3/4**, which agree among themselves to 0.109 px. A pure translation: a
rotation/drift term is not supported, since the translation-only model already
sits at the 0.07 px noise floor. 0.8115 px is the same quantity as the 0.7938 px
measured against the official i2d - the two numbers are the same error seen
twice. (Those are diagnostic values; the corrections actually *shipped* are
re-gauged to visit 1 and are listed under "Measured Stage 4 run" above. The
0.8115 px separation and the 0.109 px spread are gauge-invariant.)

**OPEN, UNTESTED LIMITATION.** The correction is applied per *visit*, which
assumes pointing is constant across all 8 detectors within a visit. That
assumption is **untested and not resolvable with the present overlap graph**,
because there are no same-visit cross-detector star matches to measure it
with. A per-detector pointing term would be invisible to the current solve and
would survive the correction. Do not describe this as resolved.

**F090W later MEASURED the size of this term** - see "F090W: the per-detector
term is real" below. It is ~0.03 px for F200W but up to ~0.27 px for F090W, and
it is the dominant residual there. The assumption is still not *corrected*,
only now *quantified*.

## F090W: the per-detector term is real

F090W (160 frames, all 8 detectors x 4 visits) was run through the same
five-step procedure. It repeated the F200W gauge answer - **visit 1** - but
that was measured, not assumed, and the two filters behave quite differently
downstream.

| filter | v1 | v2 | v3 | v4 | within-visit internal rms per visit |
|--------|----|----|----|----|--------------------------------------|
| F200W | 0.004 | 0.721 | 0.844 | 0.940 | 0.031 / 0.032 / 0.037 / 0.031 px |
| F090W | **0.017** | 0.528 | 0.808 | 0.925 | 0.056 / **0.274** / 0.167 / 0.057 px |

(Raw offset of each visit's header WCS versus the official i2d, in i2d px,
median over the 8 detectors. The F200W row reproduces the table above, which
validates the method.)

**F200W's 8 detectors agree with each other to ~0.03 px within a visit.
F090W's do not** - up to 0.27 px rms in visit 2. Because the detectors tile
disjoint sky they cannot be compared to each other, but each can be compared
to the i2d, which measures the same thing. In visit 2 the split is stark:
nrca1-4 sit at dy ~ +0.67, nrcb3/b4 at ~ +0.59, but nrcb1/b2 at ~ +0.13 - a
~0.5 px per-detector step.

This is **not** a positional gradient. Regressing each detector's offset
against its field centre gives slopes of 1e-5 to 2e-5 px/px, i.e. only
0.08-0.16 px across the 7800 px field, which cannot account for 0.27 px of
scatter. Nor is it uniform across visits (v1 and v4 are tight at 0.056 px, v2
and v3 are not), so it is not a fixed per-detector bias either. It behaves
like visit-dependent, detector-level pointing error.

### What that does to the result

| quantity | F200W | F090W |
|----------|-------|-------|
| median dx, dy vs i2d | -0.0240, +0.0327 | -0.0195, +0.0049 |
| median \|offset\| | 0.0730 | **0.1384** |
| offset MAD | 0.0430 | 0.1170 |
| offset 16-84 pct | - | 0.061 / 0.652 |
| within 0.5 px | 99.0% | **78.5%** |
| matched stars | 5882 | 5428 |
| source-free RMS (MJy/sr) | 0.0873 | 0.0609 |
| background difference (MJy/sr) | +0.0431 | +0.0286 |

Both filters end with **zero global bias** - the gauge is right in both, and
F090W's median dx/dy is if anything better than F200W's. The difference is
entirely in the *spread*: F090W's 16-84 percentile range of 0.061/0.652 px is a
tight core with a long tail, which is the signature of per-visit, per-detector
error rather than a single global shift. 78.5% within 0.5 px is the cost.

The mechanism is confirmed by how far the fitted correction departs from the
anchor-referenced separation measured in step 2:

| visit | step 2 separation from v1 (i2d px) | fitted \|shift\| (grid px) | difference | internal rms |
|-------|-------------------------------------|--------------------------|------------|--------------|
| 2 | 0.5251 | 0.7512 | **+0.2261** | 0.274 |
| 3 | 0.8046 | 0.8416 | +0.0370 | 0.167 |
| 4 | 0.9105 | 0.9240 | +0.0135 | 0.057 |

The mismatch tracks the visit's own internal scatter almost exactly. Where
pointing is rigid (v4) the rigid fit is right; where it is not (v2) the fit
overshoots, because no single translation can represent a distorted visit. This
is the open limitation above, showing up in the final residual.

**Do not report the F090W 0.1384 px as a failure of the solve.** The
pre-registered prediction was 0.1-0.2 px and the measured value is 0.1384 px;
the anchor sits 0.017 px from the i2d, exactly as step 2 predicted, and the
median dx/dy is zero. What the number measures is a per-detector pointing term
that the current per-visit model does not (and, with no same-visit
cross-detector star matches, cannot) remove. Fixing it needs a per-detector
translation term and a way to solve it - which requires either same-visit
cross-detector matches or an external reference that already encodes
per-detector geometry.

> **0.1384 px is a real, currently-irreducible limit for a per-visit
> translation model on this filter's detector geometry - not a solve failure,
> and not fixed by correcting to the i2d, which would be tautological.**

This is the shipped, intended F090W result: the independent reduction, gauged to
an externally measured anchor, validated against the official i2d. Do not
"improve" it by importing per-detector offsets from that same i2d. The
overlap graph is structurally incapable of identifying a per-detector term
anyway: every edge joins two *different* visits **and** two *different*
detectors (0 same-visit edges, 0 same-detector edges; rank 17 of 32 nodes), so
visit and detector offsets are confounded and 15 of 32 per-group parameters are
unconstrained. Measured per-group precision is good - 0.009 px median error from
65-498 matched stars per group - so the obstacle is identifiability and
independence, not noise. Correcting to the i2d would drive the reported offset
to ~0 *by construction* and destroy the only evidence that this pointing error
exists, which is the same mistake as gauging to the consensus.

### F090W measured run

```
.venv\Scripts\python -m jwst_stack.cli mosaic `
  --detector nrca1 nrca2 nrca3 nrca4 nrcb1 nrcb2 nrcb3 nrcb4 --filter F090W `
  --grid out\grid.fits --out out\f090w_all_detectors_mosaic.fits `
  --registration out\mosaic_registration_f090w.json --scratch-dir out `
  --recompute-registration --gauge-visit 1
```

| quantity | value |
|----------|-------|
| frames combined | 160/160 |
| grid | shared ny 15895 x nx 22130 (unchanged) |
| tiles | 352 at 1024 px |
| covered pixels | 120,288,441 (34.2% of grid) |
| max depth | 17 |
| per-group registration | median \|shift\| 0.048-0.376 px by group |
| cross-visit stage | 23 edges, residual rms 0.0895 px, gauge v1 (explicit) |
| corrections | v1 (0, 0) v2 (+0.5012, +0.5596) v3 (+0.4404, +0.7172) v4 (+0.5872, +0.7134) |
| frames shifted | 128/160 |
| wall time | ~31 min build |
| peak RSS | 4.58 GB (compare: 5.40 GB) |
| validation wall time | 4.0 min |

Note the F090W group shifts are larger than F200W's median 0.1284 px (v1 groups
alone sit at 0.30-0.38 px), which is the same per-detector effect showing up
before the cross-visit stage.

### Where the 0.1384 px actually comes from

Split by fitting each group's matched stars against their own i2d position
(`out/f090w_step2_detectors.json`, `out/f090w_within_detector.json`). The fit
must be sigma-clipped: per-star RMS is ~0.7 px against a MAD of 0.11 px, so plain
least squares is set by a handful of blended and saturated cores.

| component | what it is | median | detectable? |
|-----------|------------|--------|-------------|
| inter-detector | the 8 detectors of a visit disagreeing on a constant offset | 0.1118 px rms | **no** - rank 17 of 32, confounded with visit |
| intra-detector | offset varying *within* one detector's field (rotation/skew, median 0.0071 deg) | 0.0927 px rms, 0.0812 px median | **yes** - a smooth function of sky position |
| the two in quadrature | | **0.1452 px** | |
| measured mosaic median \|offset\| | | **0.1384 px** | |

The two components bracket the shipped number closely - 0.1452 against 0.1384,
5% high - so the residual is *mostly* accounted for. Do not quote this as an
exact closure: it is a consistency check on two independently measured terms,
not a derivation, and F200W below shows it is not exact. Per visit the
intra/inter ratio is 0.80x (v1), 0.49x (v2), 1.01x (v3), 0.79x (v4). The field
term is significant at a median **6.7 sigma**, with 29 of 32 groups above
3 sigma, so it is not noise.

**The same measurement on F200W** (`out/f200w_within_detector.json`): inter
0.0337 px, intra 0.0401 px rms (0.0332 median), quadrature 0.0524 px against a
measured 0.0730 px. So **F200W is not clean** - it carries an intra-detector
field term of 0.0332 px at a median **9.9 sigma** (31 of 32 groups above
3 sigma), essentially equal to its inter-detector term. The two do not cancel;
they add. The two components account for roughly 70% of F200W's residual and
leave about **0.05 px unexplained**, which is why the F090W agreement to a few
per cent must not be read as a complete budget. Candidates for the leftover
include the per-star centroid floor and PSF-difference biases against the
drizzle; it is not characterised.

Note the practical consequence: F200W and F090W fail differently. F090W is
dominated by a large inter-detector term (0.1118) that no internal measurement
can constrain; F200W's terms are each about 0.03-0.04, smaller but not zero.

**This is not an i2d rotation artefact.** A global field rotation in the i2d
would give the same trend in all 32 groups. Instead the trend ranges over
0.011-0.240 px and tracks the visit (v1/v4 tight at ~0.045 px, v2/v3 large at
0.13-0.17 px), i.e. it is per-group and per-visit, like the inter-detector term.

**So a genuinely non-i2d-dependent fix does exist, and it is worth about 20%.**
Unlike a per-detector step, a field gradient is identifiable from internal data
alone: the 5 dithers of each group overlap the same sky, so each frame's own
stars supply the positional lever arm. The fix is to let the per-group
registration fit a rotation/skew as well as a translation, which the current
translation-only model cannot represent. Stage A's own residual rms (0.125-0.438
px per group) is the same order as the corner-max field term (0.03-0.75 px),
which corroborates that the translation-only dither fit is absorbing a rotation
it cannot express.

Be precise about the payoff, though: it removes the 0.0812 px term and leaves
the 0.1118 px inter-detector one, so F090W would go from 0.1384 to about
**0.11 px**, not to zero. The larger half is the unidentifiable part.

## F187N: the third 160-frame filter, and the best result

F187N is the reproduction case the gauge procedure was written for, and it
behaves better than either earlier filter. Step 2 was done properly - all 8
detectors per visit, threshold 20 matches - before solving, and **no
low-match fallback was needed**: all 32 groups cleared the threshold with
109-604 matched stars, so every per-visit median is an 8-detector median.

| visit | dx (i2d px) | dy | \|offset\| | spread | internal rms | n_det |
|-------|-------------|----|--------|--------|--------------|-------|
| v1 | +0.0003 | +0.0060 | **0.0060** | 0.0247 | 0.0199 | 8 |
| v2 | -0.3471 | +0.5682 | 0.6658 | 0.0518 | 0.0348 | 8 |
| v3 | -0.5274 | +0.6059 | 0.8033 | 0.0627 | 0.0428 | 8 |
| v4 | -0.6335 | +0.7248 | 0.9626 | 0.0596 | 0.0254 | 8 |

**Anchor: visit 1, measured at 0.0060 px.** Note this coincides with the
min-visit heuristic, which is exactly the situation the heuristic cannot be
trusted in - it is right here because step 2 measured it, not because visit 1
happens to be lowest. The inter-visit pointing is the largest of the three
filters (v4 at 0.9626 px), and it is the *cleanest*: within-visit internal rms
is 0.020-0.043 px, versus F090W's 0.056-0.274. F187N's per-visit pointing is
effectively rigid, which is why the fit tracks the anchor-referenced separation
closely (differences of only +0.043, +0.053, +0.011 px for v2/v3/v4, against
F090W's +0.226 for v2).

### F187N measured run

```powershell
.venv\Scripts\python -m jwst_stack.cli mosaic `
  --detector nrca1 nrca2 nrca3 nrca4 nrcb1 nrcb2 nrcb3 nrcb4 --filter F187N `
  --grid out\grid.fits --out out\f187n_all_detectors_mosaic.fits `
  --registration out\mosaic_registration_f187n.json --scratch-dir out `
  --recompute-registration --gauge-visit 1
```

| quantity | value |
|----------|-------|
| frames combined | 160/160 |
| grid | shared, reused unchanged (see the axis-order gotcha below) |
| tiles | 352 at 1024 px (max 27 frames/tile) |
| covered pixels | 120,286,701 (34.2% of grid) |
| max depth | 17 |
| per-group registration | median \|shift\| 0.1030 px, 128/160 frames shifted |
| cross-visit stage | 24 edges, residual rms **0.0248 px**, gauge v1 (explicit) |
| corrections | v1 (0,0) v2 (+0.5030,+0.4996) v3 (+0.4723,+0.7140) v4 (+0.5763,+0.7846) px |
| wall time | 34.2 min build |
| peak RSS | 4.47 GB |
| validation wall time | 4.1 min, peak RSS 5.41 GB |

### F187N validation vs the official combined i2d

| quantity | value |
|----------|-------|
| overlap | 120,281,685 px (34.2% of mosaic, **97.7% of i2d**) |
| background (median) | mosaic 15.6488 vs official 14.9593 -> **+0.7185 MJy/sr** |
| residual noise, whole overlap | MAD 0.3494, RMS 32.2285 MJy/sr |
| residual noise, source-free (96.9% of px) | median 0.7187, MAD 0.3438, RMS 0.8104 MJy/sr |
| stars mine/i2d | 6000 / 6000 (both hit the `max_stars` cap) |
| matched stars | 5908 |
| star offset (mine - official) | median dx **-0.0070**, dy **+0.0215** px; median \|offset\| **0.0647 px** |
| offset scatter | MAD 0.0378 px, rms dx 0.1546 / dy 0.2510 px, 16-84 pct 0.032/0.110 |
| within 0.5 px | **99.3%** |

The prediction was registered **before** building: "0.05-0.10 px, centred ~0.07",
from F187N's measured inter-detector term of 0.0290 px being comparable to
F200W's 0.0337 rather than F090W's 0.1118. Measured 0.0647 px, inside the range
and near the centre. The gauge is right (median dx/dy within 0.02 px of zero).

Two F187N-specific numbers that are **not** astrometric and should not be read
as regressions:

- The **background difference is +0.7185 MJy/sr**, ~20x F200W's +0.0431 and
  F090W's +0.0286. F187N is a narrowband filter over a bright PAH continuum
  (14.96 MJy/sr official background, 5.4x F200W's), so this is a photometric
  calibration difference of the same kind as the other two, just scaled by the
  background. As with F200W, a rigid star translation cannot move a median
  background, so this is not an alignment artefact.
- The **source-free RMS (0.8104) and MAD (0.3438) are correspondingly high**,
  and are dominated by that same global offset: the source-free *median*
  residual is 0.7187, i.e. essentially the entire offset. F187N's noise floor
  is not comparable to F200W's 0.0873 in absolute terms.

### The three filters together, and a ~0.05 px common floor

| filter | inter-detector | intra (median / RMS) | in quadrature | measured median \|offset\| | variance explained | field-term significance |
|--------|---------------|----------------------|---------------|---------------------------|---------------------|------------------------|
| F090W | 0.1118 | 0.0812 / 0.0927 | 0.1452 | **0.1384** | 110% (over) | 6.7σ, 29/32 |
| F200W | 0.0337 | 0.0332 / 0.0401 | 0.0524 | **0.0730** | 52% | 9.9σ, 31/32 |
| F187N | 0.0290 | 0.0263 / 0.0316 | 0.0429 | **0.0647** | 44% | 6.0σ, 29/32 |

Read this carefully, because it is easy to over-claim:

- F187N and F200W are the same shape of problem. Both carry a real intra-detector
  field term of the same order as their inter-detector term (ratios 0.91x and
  0.98x pooled, 0.84-0.96x per visit for F187N), and both **fail to close**:
  subtracting the quadrature leaves **0.0485 px (F187N) and 0.0508 px
  (F200W) unexplained**. Those two leftovers agree to 5% despite the filters
  differing in brightness by 5x, which suggests a common systematic of about
  0.05 px in this comparison rather than noise - candidates remain the
  per-star centroid floor and PSF-difference bias against the drizzle. It is
  **not characterised**, and it is a residual, not a measured physical term.
- F090W does **not** show this floor, and cannot: its measured 0.1384 is
  already *below* its 0.1452 quadrature, so there is nothing to subtract. Its
  110% "variance explained" means the two terms slightly over-predict, which is
  the consistency check succeeding, not a physical claim.
- So the honest summary is a floor, not a ranking: **no per-visit translation
  model gets below ~0.05 px against these i2d products**, and F090W's 0.1384
  is dominated by a separate, larger, unidentifiable per-detector step.
- The field term is significant in all three filters (29-31 of 32 groups >=3σ),
  so it is not noise anywhere. It is real, it is measurable, and it is the
  part the deferred rotation-aware fit would remove.

### The gauge is the part that is easy to get wrong

The correction is a *relative* translation, so its absolute origin is
unobservable from the star matches alone and must be chosen. The first
implementation gauged to the **consensus** of mutually-agreeing visits
(v2/3/4), and it was wrong. It produced a mosaic that was internally consistent
(offset MAD collapsed 0.2275 -> 0.0613 px, proving the inter-visit error was
genuinely removed) yet sat 0.86 px from the official i2d - exactly the size of
the correction just applied.

Why: visits 2/3/4 agreeing with each other says nothing about where the mosaic
sits in absolute terms. Gauging to them silently redefines the sky frame.

Measured against the official i2d, the **raw, uncorrected** frame WCSs sit:

| visit | offset vs i2d (i2d pixels) | \|offset\| |
|-------|---------------------------|-----------|
| 1 | +0.0097, -0.0031 | **0.0101 px** |
| 2 | +0.4087, -0.5926 | 0.7199 px |
| 3 | +0.6178, -0.5839 | 0.8501 px |
| 4 | +0.6137, -0.7048 | 0.9346 px |

So the official i2d agrees with **visit 1's header WCS to 0.010 px**. The
gauge is therefore pinned to the lowest-numbered visit, which leaves that visit
exactly where the header WCS put it and moves the others onto it. This also
independently confirms the solve: rotating the internal grid-pixel corrections
by the -102.43 deg grid-to-i2d rotation reproduces these i2d-frame numbers to
~0.02 px on every visit.

**Never gauge a relative astrometric correction to a consensus cluster**
without checking what the consensus means in absolute terms. The consensus is
still computed, but only as a diagnostic for spotting a displaced visit.

### The rule that generalises, and the per-filter procedure

> Internal self-consistency among a majority does not imply external
> correctness. Always cross-check the gauge choice against an independent
> reference before trusting it.

This is not an F200W quirk. It applies to **every** filter, and the trap
recurs with **different visits playing the "consensus" role**. For F200W the
displaced visit happened to be visit 1 and the anchor was also visit 1, so the
min-visit heuristic gave the right answer *by luck*. A filter where visits
1/2/3 agree and visit 4 is displaced would have the min-visit heuristic pin
the **wrong** frame - and the run would look healthy: same 0.07 px internal
residual, same collapsed offset MAD, same green validation, sitting a full
inter-visit error away from the reference.

**The anchor is a per-filter measurement, not a convention.** Nothing in the
data identifies it, and no internal consistency check can.

Procedure, per filter, before trusting a mosaic:

1. Check whether the existing fixed grid already contains the filter's padded
   footprint. `out/grid.fits` is a **program-wide** grid, not per-filter: the
   union of *all* 160 F090W footprints fits inside the F200W grid, as does
   F187N's (x 5303..16806, y 24..15873, margins 5303/5323/24/21 px), so the
   grid is reused as-is. Do not rebuild it per filter. Do the check with
   `_grid_bbox(exp, grid_wcs, shape)`, which returns `(y0, y1, x0, x1)` and
   takes the *exposure*, not a WCS. The **grid-to-i2d similarity transform and
   the gauge are per filter** and must never be reused from another filter.
2. Star-detect one reference frame per visit, convert to i2d pixels, and
   cross-match each visit against a common reference. This measures the **raw,
   uncorrected** frame WCS offsets versus the official i2d, per visit. Do this
   for **all 8 detectors** per visit, not one: a single reference frame is
   several times noisier (F090W: 0.119 px from one nrca1 frame versus 0.017 px
   from the 8-detector median), and the per-detector spread is itself the
   quantity of interest (see the F090W section below).
3. The visit whose offset is ~0 is the anchor. That is `--gauge-visit`.
4. Only then solve, passing `--gauge-visit` explicitly, and rebuild.
5. Validate against the i2d as normal. The expected median |offset| is then the
   *predictive* residual of step 2, not 0. **State the predicted range before
   building.** It is not always ~0.01 px: it is the per-visit scatter that a
   rigid translation cannot remove, so it differs sharply per filter (F200W
   0.073 px measured, ~0.01-0.02 px predicted; F090W 0.138 px measured, 0.1-0.2
   px predicted). If the measured value does not match the prediction, the
   solve is wrong - do not explain it away.

`gauge_source` in the registration JSON records `explicit` vs
`heuristic:min-visit`, and the solve prints a warning when a displaced visit is
present but the gauge was not pinned externally. Treat
`heuristic:min-visit` as unvalidated until step 2 is done. The shipped F200W
`out/mosaic_registration.json`, F090W `out/mosaic_registration_f090w.json` and
F187N `out/mosaic_registration_f187n.json` all record `gauge_source: explicit`
with `gauge_visit: 1`, because step 2 was done for each of those filters
(0.0101, 0.017 and 0.0060 px).

An external reference exists for **every** NIRCam filter, so step 2 is always
possible. Combined all-visit products in program 2731:

| filter | combined i2d | cal files | full 8-det x 4-visit? |
|--------|---------------|-----------|----------------------|
| F090W | `jw02731-o001_t017_nircam_clear-f090w_i2d.fits` (5,415,445,440 B) | 160 | yes - **done** |
| F187N | `jw02731-o001_t017_nircam_clear-f187n_i2d.fits` (5,416,456,320 B) | 160 | yes - **done** |
| F200W | `jw02731-o001_t017_nircam_clear-f200w_i2d.fits` (5,420,381,760 B) | 160 | yes - **done** |
| F335M | `jw02731-o001_t017_nircam_clear-f335m_i2d.fits` (944,968,320 B) | 40 | no, single visit |
| F444W | `jw02731-o001_t017_nircam_clear-f444w_i2d.fits` (944,968,320 B) | 40 | no, single visit |
| F444W;F470N | `jw02731-o001_t017_nircam_f444w-f470n_i2d.fits` (944,254,080 B) | 40 | no, single visit |

All three 160-frame filters are downloaded and done. The F335M/F444W sets are
single-visit and have no inter-visit alignment to solve at all. (MIRI
F1130W/F1280W/F1800W/F770W are a different instrument and out of scope for the
NIRCam code path.)

### Gotchas when re-measuring astrometry here

- **`load_grid` returns the shape as `(ny, nx)`, and the grid file has no
  `NAXIS1`/`NAXIS2`.** `load_grid` reads `GRDNX`/`GRDNY` and returns
  `(nax2, nax1)`, so `shape[0]` is the *height*. The grid is ny 15895 x
  nx 22130, i.e. **wider than tall**; `out/f200w_all_detectors_mosaic.fits` has
  `data.shape == (15895, 22130)`. Reading `shape` as `(nx, ny)` makes a valid
  footprint union look like it overruns the grid: F200W, the filter the grid
  was actually built from, "needed" 16812 x 15874 against a "grid" of
  15895 x 22130 and appeared not to fit. The real footprint is x 5316..16796,
  y 36..15858 (read from the finished mosaic's nonzero SCI), with ~5300 px
  margins in x and 36 px in y. Unpacking it backwards would have rebuilt the
  grid and invalidated all three shipped mosaics. Probing `header["NAXIS1"]`
  raises `KeyError` for the same reason - the keywords are `GRDNX`/`GRDNY`.
- **A WCS-only per-frame residual test is tautological.** Mapping
  `frame pixel -> sky -> grid` and `sky -> i2d`, then comparing, returns
  exactly 0.0000 px for any frame WCS whatsoever, because the frame WCS is
  only a coordinate generator and cancels. Measured over 46,240 points from
  all 160 frames: 0.0000. It cannot detect frame-level astrometric error; that
  requires stars. Do not use it as evidence.
- **The grid-to-i2d relation is an exact similarity transform**: scale
  0.992731 (= our 0.0310"/px / i2d 0.031227"/px), rotation -102.4292 deg, max
  fit residual **9.24e-05 px** over 682 points, unit anisotropy. Any apparent
  "TAN vs drizzle convention offset" must be checked against this first - the
  reprojection already handles the rotation and scale exactly, so a residual
  offset is ours, not the convention's. **Vectors measured in grid pixels are
  not in i2d pixels**: rotate by -102.43 deg before comparing, or the signs
  will look wrong.
- **Header CRVALs cannot corroborate inter-visit offsets.** The four visits
  point 40-260" apart and share no like-for-like dither cells, so mean CRVAL
  per visit mixes pointing with dither geometry. Star matching is the only
  route.
