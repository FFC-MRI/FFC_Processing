# Reconstruction loading hotfix — 7 September 2026

This update only changes Cameleon stored-average handling, preprocessing
diagnostics and failed-load handling. The existing fitting, SVS reduction,
orientation and voxel-localisation implementations are unchanged.

## Cause and correction

The Cameleon reader returns `[read phase slice experiment receiver]`. Its
`par.navg` describes scanner accumulations (`NUMBER_OF_AVERAGES`), not an
additional stored dimension. The constructor had incorrectly copied that
count into the dimension used by preprocessing.

For the supplied sequence 1, `128 * 128 * 15 * 13 = 3,194,880` complex samples
are stored (25,559,040 bytes of float32 IQ data). `NUMBER_OF_AVERAGES=6` made
preprocessing request six times that number. The `.dat` constructor now sets
the stored average dimension to one and preserves the acquisition count in
`scan.param.NUMBER_OF_AVERAGES`. It does not divide the signal by six. Legacy
formats with explicit stored average blocks still use the existing mean over
that dimension.

The later receiver error was secondary: the GUI retained a scan after its
preprocessing failed, then passed empty k-space to reconstruction. MATLAB
reports `size([],6)` as one, hence the misleading "available coils 1..1".

## Loading behaviour

- Imaging preprocessing validates and reports both actual sample count and
  expected dimension sizes before reshaping; it does not infer a replacement
  coil count or truncate data to conceal a mismatch.
- Reconstruction reports `ImageReconCore:NoPreprocessedKspace` if it is called
  without successfully preprocessed imaging data.
- The loader keeps only successfully preprocessed scans, with matching names
  and selection indices. Partial failures are reported with full stack traces
  in the Command Window and a warning in the app.
- Cancellation or an entirely failed batch leaves an already loaded dataset
  unchanged and does not call reconstruction on the failed input.

## Installing and checking

Close FFCProcessing and restart MATLAB so that the old class/app code is no
longer cached. Extract this ZIP into the existing FFC_Processing folder,
replacing matching files and retaining all subfolders. Reopen FFCProcessing
and reload the original acquisition files (not a cached processed object
created by the faulty version). No header or raw data edits are needed.

The changed runtime files are `ImageReconCore.m`, `FFCProcessing.mlapp`, and
the new `ffc_preprocess_scans.m`. Install them together.

From the project root, run:

```matlab
results = runtests('tests');
table(results)
assertSuccess(results)
```

`tests/test_cameleon_loading.m` covers 2D/3D data, one/thirteen receivers,
one/six scanner accumulations, the reported 128-by-128-by-15 matrix, image
reconstruction, time/field ordering, SVS magnitude reduction, explicit stored
average blocks, invalid dimensions, absent preprocessing, real invalid coil
selections, mixed-success batches, all-failed batches and empty input.

## Validation limits

MATLAB is not installed in the editing environment, so the MATLAB tests and
the App Designer runtime have not been executed here. Independent numerical
checks passed for all five supplied data/header pairs, including the new
sequence 1: stored sample counts, dimension mapping and every receiver's
complex samples were preserved exactly. The former sixfold reshape failure
was also reproduced on sequence 1.

Static MATLAB syntax checks passed for the modified class, batch helper,
new tests and app callback source. The packaged app's XML was checked, and
unchanged archive members (including all fitting and geometry helpers) were
compared byte-for-byte with the previous delivered ZIP. MATLAB-side tests
are included for verification on the acquisition workstation.
