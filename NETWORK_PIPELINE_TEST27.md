# TEST27 — network-drive loading pipeline

Built from TEST26 on 10 November 2026. Numerical reconstruction, geometry,
DICOM export and TEST26 cache validity rules are unchanged.

## Start

Close the old FFC app, extract to a separate folder, set MATLAB's current folder
to FFC_Processing_TEST27, and run `ffc_start_cvi`.

Folder loads and raw-only DAT selections now use one client reader and up to two
concurrent reconstruction jobs. The client reads each sequence and sends the
loaded object to a process worker for preprocessing and reconstruction. The next
network read can overlap this work. Completed results are collected in original
sequence order; every successful sequence is ready before the viewer opens.

At most two jobs are in flight by default. There is no unbounded queue of raw
sequences. Completed images and their source data still remain in RAM as in
TEST26; this is not disk paging or a total-study RAM limit. Process workers
require additional memory for transferred data and reconstruction temporaries.

Workers perform no raw-data or cache IO. Cache lookups/writes remain serial on
the client. A valid local cache hit bypasses worker transfer and pool startup.
The existing TEST26 cache is deliberately retained and remains compatible.

A process pool is started only on the first uncached scan of a multi-sequence
load. Existing process pools are reused with at most the configured number of
jobs; the pool is retained for future loads. Pool startup can add a noticeable
one-off delay. Without Parallel Computing Toolbox, on pool-startup failure, or
with an existing thread pool, loading falls back to serial reconstruction and
reports why. The user's existing pool is not closed or replaced. Single-sequence
loads run serially because there is no following read to overlap.

MAT/MRD and mixed-format selections keep the TEST26 serial loading path.

## Progress, failure and cancellation

The progress dialog shows current stage, sequences processed (ready or failed),
active jobs, elapsed time, raw/header reading time, and completed reconstruction
work. Failures are reported per sequence; other sequences continue. Cancelling
cancels outstanding futures and preserves the previously displayed image set.
A synchronous network read, pool startup, or local cache write must return before
the client can observe cancellation. Serial fallback cannot interrupt an FFT.

The MATLAB Command Window prints elapsed/read/cache/preprocess/reconstruction
summaries. Worker times are sums of completed job times and overlap: they must
not be added together to infer total elapsed time. Wall time also includes
worker transfer and polling overhead. Folder discovery/file-picker time is not
included. Raw/header timing combines source validation, reading, parsing and
construction; this version does not separately time XML parsing versus binary IO.

For detailed results, launch with an output argument:

```matlab
app = ffc_start_cvi;
% Open a raw image set using the app, then:
report = app.getLoadTimingReport();
disp(struct2table(report.sequence))
```

## Controls

Defaults are parallel enabled, two jobs. Changes apply to the next load:

```matlab
setpref('FFCProcessing','ParallelLoading',true)
setpref('FFCProcessing','LoadWorkers',2)
% To compare against serial reading/reconstruction:
setpref('FFCProcessing','ParallelLoading',false)
```

Use the same study and Reprocessing setting when comparing. Reprocessing ON
bypasses the disk cache; OFF tests normal cached reopening. Restore parallel
loading afterwards. Increasing worker count may increase RAM use without helping
when the network is already the limiting factor. Existing larger pools are not
resized; only the number of submitted jobs is capped.

## Validation

Changed MATLAB files pass MISS_HIT parsing/static checks. MATLAB is unavailable
in the build environment: runtime correctness, UI responsiveness and real-world
speedup have NOT been measured. New tests in tests/test_network_pipeline.m cover
serial fallback, failure isolation/order, cancellation, warm-cache worker bypass,
parallel versus serial image equality, and the two-job bound. The actual parallel
tests require Parallel Computing Toolbox and a process pool (or permission for
MATLAB to start its local pool).

```matlab
ffc_setup
results = runtests('tests');
disp(table(results)); assertSuccess(results)
```

Suggested real-data check: load a network study with Reprocessing ON; verify
all sequences open ready, ordering/geometry match TEST26, and the progress report
shows active reconstruction jobs while subsequent files are read. Repeat with
parallel loading disabled for timing comparison, then reopen with Reprocessing
OFF to check the local cache. Cancel a load mid-way and verify the previous study
remains visible. Full TEST26 behaviour and cache preferences are documented in
PERFORMANCE_TEST26.md (historical baseline).
