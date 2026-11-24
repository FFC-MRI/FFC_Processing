# TEST28 — all loaded participants and parallel DICOM export

Start `ffc_start_cvi` from the extracted FFC_Processing_TEST28 folder after
closing the older app. The window title identifies TEST28.

The Enhanced DICOM export button now offers:
- Selected sequence (the existing single-series workflow).
- All loaded participants (replaces the selected-participant batch option).
- Cancel.

Load the common parent folder containing the participants' raw datasets first.
All imaging sequences currently in the scan list are then included in the batch;
the highlighted participant does not restrict the export. Opening a new dataset
still replaces the loaded image set rather than appending to it.

Select a destination folder. Each export creates a fresh FCI_export_<UID>
subfolder containing all participants' DICOMs together. Import that subfolder
into CVI. Filenames include scan index, patient label and series description;
repeated series names/numbers cannot overwrite another participant's images.

Identity and grouping:
- Metadata is extracted independently for each scan on the MATLAB client.
- Each worker receives that metadata snapshot, avoiding repeated network XML
  reads and any borrowing from the selected participant.
- Missing name, DOB or sex can be filled only from a unanimous value among scans
  with the same explicit PatientID. Existing non-empty values are never replaced.
- Patient ID scopes study/frame-of-reference identifiers; name and DOB are the
  fallback when ID is absent. Study ID, description and date separate visits.
- If the study date is missing, studies are kept separate per scan instead of
  guessing that they are one visit. Batch study UID generation now includes
  explicit patient identity and differs from the historical single-export scheme.
- Scans lacking both patient ID and name are skipped with a diagnostic; identity
  is not invented. SVS/non-image sequences are also skipped.
- The DICOM geometry, intensity scaling and frame order are unchanged.

Parallel export defaults to two process-worker jobs. An existing process pool
is reused; otherwise a local pool is started. Without a suitable pool/toolbox,
export runs serially. Valid image snapshots omit raw k-space and unused image
arrays before transfer to reduce worker memory/transfer cost. Stale images are
rebuilt as needed and returned to the viewer. There is no nested parallel work.
Each file is written to a temporary .part file and published only when complete.

Progress reports completion/failures. Cancellation stops outstanding jobs and
keeps completed files. A job that finishes as cancellation is requested may have
published a file before its result is collected; the cancelled summary counts
collected results. Pool startup/serial writing must return before cancellation
can be observed. Failed/skipped series are listed in the Command Window.

Optional MATLAB preferences (apply to the next export):

```matlab
setpref('FFCProcessing','ParallelDicomExport',true) % default
setpref('FFCProcessing','ExportWorkers',2)          % default
% Compare serial export, particularly for a slow network destination:
setpref('FFCProcessing','ParallelDicomExport',false)
```

Parallel writing may be limited by destination/network throughput. No measured
speedup is claimed. TEST27 network loading and the TEST26 local cache are retained.

Validation: changed files pass MISS_HIT parsing/static checks. MATLAB is not
available in the build environment. Runtime tests and CVI import validation
remain outstanding. Tests cover cross-participant/visit UID separation, unique
filenames, same-ID demographic completion, anonymous/SVS exclusions, DICOM
metadata readback, and parallel versus serial pixel equivalence.

```matlab
ffc_setup
results = runtests('tests');
disp(table(results)); assertSuccess(results)
```

The readback/parallel tests require the applicable MATLAB toolboxes. For a real
smoke test, load two participants with duplicate sequence names, export all,
and confirm CVI shows the correct patient ID, DOB/sex and study for each.
