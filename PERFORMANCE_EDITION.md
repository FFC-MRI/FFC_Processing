# FFC Processing — Performance Edition

Prepared 7 September 2026 from `FFC_Processing-master(8).zip`.

This version implements reconstruction reuse, bounded reconstruction intermediates, lighter display updates, consistent cached filter previews, and managed background jobs. MATLAB runtime validation and timings are still required on your machine. The checks completed here are recorded below.

## Start the app

Extract the ZIP into a separate folder. In MATLAB, make `FFC_Processing` the current folder and run:

```matlab
ffc_start
```

The active app is now the editable `FFCProcessing.m` class. It contains the existing interface and the updated callbacks. The original App Designer file is preserved in `legacy_app/FFCProcessing_original.mlapp` for reference. The new source does not depend on editing or synchronising App Designer's binary component model.

Use `ffc_setup` to configure the active paths. Avoid combining this folder with another FFC installation on the MATLAB path. Restart MATLAB when switching between the old and new class definitions.

The target is MATLAB R2025b, matching the uploaded app's metadata. Keep the toolboxes required by your existing workflows; Wavelet filtering also requires Wavelet Toolbox. Background work uses `backgroundPool`. No process pool is started when loading files or opening a filter preview.

## What changed

| Area | New behavior |
|---|---|
| Scan selection | Reuses a valid reconstruction and preserves filtered/registered magnitude data, slice, time, field and window settings. |
| Reconstruction | Whitens and reconstructs one time/field volume at a time. All partitions of a true 3D volume remain together. |
| Calibration | Retains the small noise covariance/factor across compatible FFT/window/coil changes; forced rebuilds recalibrate. |
| Zero filling | Allocates the final centred padded array once. Zero padding and a full-size `None` window have fast paths. |
| Coil products | The viewer releases per-coil RSSQ images after combination. The public core API retains them by default for compatibility. |
| Window dragging | Changes the existing image object's `CLim`; it does not re-extract images, reconstruct, or redraw k-space. |
| Slice/time/field navigation | Updates existing image objects. K-space is cached independently; the compass redraws only when its basis changes. |
| Anatomical planes | Selects a thin native slab before permuting it, avoiding an unnecessary volume permutation. |
| Filtering | Debounces requests, caches the current preview stack, and uses identical grouping for Preview and Apply. |
| LLR | Uses one 2D coverage array, computes each patch mean once and avoids constructing a diagonal shrinkage matrix. |
| T1 maps | Converts only bounded blocks of selected voxels to double and supplies sliced inputs to optional foreground `parfor`. |
| ROI fitting | Converts selected ROI values to double rather than duplicating the whole magnitude dataset. |
| Long operations | Reconstruction, filter application, mapping and registration computation use a managed background job. |
| Loading | Removes unused pool startup and duplicate MAT loads; folder discovery reads each `data.dat` dataset once. |
| Distribution | Removes 11 autosave/log files and leaves experimental/demo folders off the default path. |

For the review example, 128 × 128 × 15 × 8 times × 9 fields × 13 receivers, a 256 × 256 padded complex-single coil array would be about 6.86 GiB if all times and fields were held together. The new padded working volume contains one time/field and is about 97.5 MiB. These are array-size calculations, **not measured peak RAM or a speedup claim**. Raw acquisition data, reconstructed products, filter output and other working arrays still consume memory.

## Processing state and numerical behavior

The source k-space remains immutable. Repeated preprocessing now reinstalls `originalcomplexkspace` instead of treating an already shifted working copy as a new source. FOV corrections remain separate from display rotations and flips. Whitening continues to use all acquired receivers before selecting virtual output channels.

The FFTs retain the first two spatial axes for 2D and the first three for 3D. Time, field and receiver dimensions are never transformed. Centred zero filling retains the existing floor-before/ceil-after convention for odd differences.

Two ACS/Walsh errors are corrected: averaging now explicitly uses dimensions 4 and 5 while keeping receivers in dimension 6, and 3D calibration images use a spatial 3D inverse FFT. **ACS/Walsh results, and ACS bias correction, can therefore differ from the uploaded version.** The covariance smoothing, principal-eigenvector calculation and phase convention remain in place. Symmetry avoids calculating both halves of the covariance independently.

The simple IR model, data-selected negative prefix, least-squares objective, optimiser bounds and tolerance are unchanged. The scalar-output fit path used for maps runs the same optimiser as the diagnostic fit. SVS still starts with mean magnitude across sampled points.

Filters are noncumulative: Apply starts from the scan's saved filter source. **Reset filter** replaces the old duplicate “Permanently apply” button and restores that source. Registration establishes a new filter source and records its masks and algorithm settings. Filtering after registration therefore uses the registered images.

Scan switching preserves the current processed magnitude. An explicit **Rebuild Images**, or a numerical reconstruction-setting change, returns to native reconstruction, reapplies the currently stored denoiser, clears registration metadata and invalidates maps/ROI fits. Rebuilding also clears the mapping mask. Changing magnitude through filtering or registration clears previous relaxation results. Maps are now empty until computed, rather than placeholder zero arrays.

The LLR grouping control offers:

- **All slices (legacy):** joins all slice/time/field frames, preserving the historical Apply behavior. The preview now uses that same full stack, so its first calculation can be expensive.
- **Per slice:** joins time/field frames within each slice. This is an explicit change to the LLR model and often reduces preview cost. Moving between times or fields of that slice reuses the result.

Stored parameters without `stack_grouping` retain legacy grouping. Wavelet, TV and TGV previews process the selected frame because those operations are independent across frames. The existing column-wise `wdenoise` behavior is retained; this release does not substitute a different wavelet algorithm.

## Background work and cancellation

The Rebuild button becomes **Cancel** during a managed task. Display navigation remains available when an existing image is present. Numerical controls are temporarily disabled, and completed results are checked against the originating dataset and magnitude generation. Switching scans, loading another dataset, closing a filter dialog or deleting the app cancels the relevant outstanding job. A failed task leaves the previously published data intact; cancelled/failed reconstruction settings return to the last successful reconstruction when one exists.

Toolbox thread support varies. If a function reports that it is unsupported in threads, the same operation is retried on the MATLAB client with a Command Window warning. That fallback is synchronous and cannot provide the same navigation/cancellation responsiveness. In particular, the existing Wavelet routine may require this fallback. General computation errors are reported, not silently converted into successful results.

Disk reads, interactive registration ROI collection, optional legacy analyses, exports, and the first automatic SVS localizer reconstruction remain client operations. A non-interruptible library call may take time to respond to cancellation. These limits are not hidden by the background-job interface.

Relevant MathWorks documentation: [backgroundPool](https://www.mathworks.com/help/matlab/ref/parallel.backgroundpool.html), [wdenoise](https://www.mathworks.com/help/wavelet/ref/wdenoise.html), [imregdemons](https://www.mathworks.com/help/images/ref/imregdemons.html).

## Validation

Completed in this environment:

- Static parsing/lint of every new or modified MATLAB source file with MISS_HIT; see `validation/matlab_static_analysis.txt`.
- 22 independent NumPy/SciPy algebra checks covering blocked whitening, centred padding, spatial FFT blocking, LLR coverage/shrinkage, plane extraction and streaming ACS means. Maximum absolute discrepancy was approximately 1.33 × 10⁻⁷. See `validation/independent_numeric_checks.json` and its reproducible script.
- Distribution integrity, source-change manifest and original-archive hash checks.

**MATLAB is unavailable here. The MATLAB tests, app launch, toolbox behavior, real-data outputs and timings have not been executed or measured.** The independent checks do not execute the MATLAB implementation.

Run the 49 MATLAB regression tests, including the existing 31 tests and 18 new numerical/state/background-lifecycle tests:

```matlab
ffc_setup;
results = runtests('tests');
disp(table(results));
assertSuccess(results);
```

One existing batch-loading test now expects empty maps, matching their deliberate lazy allocation. The new numerical tests compare blocked reconstruction with a bulk reference for 2D/3D, whitening on/off, and multiple combination modes. Separate checks exercise corrected ACS indexing, preprocessing idempotence, filter grouping, fit equivalence, and cancellation/replacement of background results.

Then check representative real scans in the app:

1. Load 2D, true 3D and SVS datasets. Check receiver IDs, FOV offsets and image orientation at each FFT size.
2. Drag window/level, scroll slices, change fields/times and return to a previously viewed scan. Check that view settings and processed images persist.
3. Apply LLR using both grouping modes. Confirm the selected preview frame matches Apply, repeated Apply is noncumulative, and Reset restores the filter source.
4. Start and cancel reconstruction/mapping/filtering; switch scans during work; close the filter dialog during a preview. Check that a late result never replaces another scan.
5. Recompute maps after filtering/registration. Compare the simple IR fit with your established results. Inspect ACS output separately because of the corrected calibration.

Measure on your own acquisition, using the same reconstruction settings for comparisons:

```matlab
ffc_setup;
scans = open_mrifile();
[scans, names, failures] = ffc_preprocess_scans(scans);
scan = scans{1};
scan.combination = 'Adaptive';
scan.noise_whitening = true;
scan.retain_coil_images = false; % matches the new viewer
report = benchmark_ffc(scan, 'ffc_timings.mat');
```

The benchmark records first reconstruction, three subsequent forced reconstructions, cached reuse and plane extraction. Compare old/new full reconstructions in separate MATLAB sessions using the same datasets and settings; cached reuse is a different operation. GUI latency and peak memory require separate measurements.

## Core API notes

```matlab
scan = buildimages(scan);        % force reconstruction (compatible default)
scan = buildimages(scan, false); % reuse if the reconstruction signature matches
scan.retain_coil_images = false; % opt out of final RSSQ coil cube
```

Scripts that directly replace raw/preprocessed arrays must call `preprocessing` with `reprocess=true`, or force `buildimages`. The reuse signature intentionally does not hash multi-gigabyte arrays. Use `ffc_invalidate_analysis` after externally modifying magnitude data. Advanced legacy fitters and the existing bias-field model remain available; coarse sensitivity maps, different solvers, reduced precision and GPU conversion were not introduced in this release.
