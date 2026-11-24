# TEST26 — reconstruct on opening and reuse

Built from the supplied TEST25(1) package on 9 November 2026.

## Start

Close the existing FFC app. Extract this package to a separate folder, make
`FFC_Processing_TEST26` the MATLAB current folder, then run:

```matlab
ffc_start_cvi
```

The title contains `Reconstruction Cache TEST 26`.

## Behaviour

- Opening a dataset preprocesses and reconstructs every successful sequence
  before publishing the new image set to the viewer. A progress dialog reports
  the sequence being processed. Cancel takes effect between sequences; a single
  reconstruction or cache write is not interruptible. Cancellation preserves
  the previously open image set. Failed sequences are reported and skipped.
- Completed images are retained in memory. Selecting another sequence only
  selects cached images; it does not run another FFT or coil combination.
- CVI masks, dispersion model fits and export metadata no longer invalidate
  reconstruction. Reconstruction keys contain explicit numerical dependencies.
- Viewer defaults are applied before reconstruction and are idempotent.
- Changing selection during an explicit rebuild/filter/map task does not cancel
  that task. Its result is committed to its originating sequence. Cancel still
  explicitly cancels the current task; rollback targets that sequence, not the
  one now selected. A stale selected scan is serviced after the current task.
- Window functions, ACS/Walsh sensitivity maps and contour intersections are
  reused when their dependencies are unchanged. Forced rebuild recalculates
  calibration/maps. Numerical coil combination and TEST25 geometry are unchanged.
- Anatomical display selects a thin slab directly from the stored 5D data.
  The viewer retains complex data for fitting and computes phase on a displayed
  plane only. Phase image exports and workspace exports materialise full phase
  arrays. Direct core callers retain the old phase/per-coil output defaults.

## Reopening raw Cameleon datasets

The raw `.dat` folder/file loader automatically attempts a local baseline cache
before reading raw data. Successful initial reconstructions are saved atomically
as versioned MAT files. Cached results include the source data needed to rebuild.
The cache is an accelerator, not a saved analysis session: later user filters,
registrations, masks and fits are not automatically saved there. Use the existing
workspace export workflow to preserve those results. MAT/MRD loads do not write
into the raw baseline cache.

Cache validation checks canonical dataset path, raw size and modification time,
exact `header.xml` and `Serie.xml` text, schema version and the reconstruction key.
A cache miss, corruption or write failure falls back to normal reconstruction.
Enable **Reprocessing** before opening to bypass the cache and reread/reprocess
raw data. Use this if raw bytes were edited while preserving BOTH file size and
modification time; this fast cache does not hash the whole raw file.

Cache location (local MATLAB preferences directory):

```matlab
fullfile(prefdir,'FFC_Reconstruction_Cache_TEST26')
```

Optional preferences:

```matlab
setpref('FFCProcessing','DiskCacheEnabled',false) % disable disk caching
setpref('FFCProcessing','DiskCacheEnabled',true)  % default
setpref('FFCProcessing','DiskCacheLimitGiB',10)   % default disk budget
% Optional alternative local cache folder:
% setpref('FFCProcessing','DiskCacheDirectory','D:\FFC_Reconstruction_Cache_TEST26')
```

Use a dedicated cache directory: eviction manages MAT files in that directory.
Oldest cache writes are removed when the disk budget is exceeded; the latest
single entry is retained even if it alone exceeds the budget. Deleting the
cache directory merely causes regeneration on the next load. Cache files contain
scan data and metadata and stay on the machine running MATLAB.

First opening now takes the reconstruction time up front and adds a cache write.
Reopening can avoid raw parsing/reconstruction, but still loads cached arrays
from disk. All reconstructed scans remain resident: this build does not evict
images from RAM or page them from disk. One-at-a-time reconstruction limits
concurrent working memory; very large image sets still require sufficient RAM.

## Validation

MATLAB is not installed in the build environment. No MATLAB runtime test, GUI
smoke test or measured speedup is claimed. The changed MATLAB files were parsed
and statically checked with MISS_HIT. See validation/TEST26_static_analysis.txt.
New regression cases are included in tests/test_performance_regressions.m for:

- Full-batch readiness and repeated A/B/A selection without new generations.
- CVI/analysis metadata preserving masks, registration and fit results.
- Real reconstruction settings/FOV offsets still invalidating cached images.
- Idempotent viewer defaults and coil-selection normalisation.
- Direct slab extraction matching previous volume extraction in all views.
- On-demand phase matching stored phase.
- Sensitivity-map reuse and receiver-change invalidation.
- Batch failure isolation and cancellation.
- Disk round trips and raw/header/Serie changes invalidating the disk cache.

Run the MATLAB suites from the new folder:

```matlab
ffc_setup
results = runtests('tests');
disp(table(results))
assertSuccess(results)
```

To time numerical batch preload and cache reuse (no GUI latency measurement):

```matlab
scans = open_mrifile();
report = benchmark_ffc_dataset(scans);
```

Suggested GUI check: load a mixed 2D/3D/SVS set, browse every sequence, perform a
CVI dispersion fit, switch away/back, and confirm the results remain. Change FFT
size and verify one rebuild; during it switch to a ready sequence and confirm
that the original completes. Check phase display/export, then reopen the raw
folder with Reprocessing off and compare with a forced reprocess.
