# TEST25: reconstruction offsets and shared patient geometry

This note supersedes the coordinate assumptions in TEST19–TEST24 and older
CVI/MPR documentation. Built from the supplied TEST24(1) source package.

## Changes

- `correct_orientation.m` subtracts the three acquisition-handled offsets in
  `OFF_CENTER_FIELD_OF_VIEW_EFF` from the prescribed local offsets before
  applying the Fourier phase ramp. Missing EFF retains legacy behaviour;
  malformed EFF is rejected. A 2-D slice index never receives a kz phase ramp.
  Offset diagnostics retain the prescribed, handled and applied values.
- `ffc_scan_spatial_geometry.m` converts Cameleon prescription vectors and
  centres with `diag([1 -1 -1])` into DICOM LPS. Cameleon is identified by the
  presence of `IMAGE_ORIENTATION_SUBJECT` in the parameter struct. Other
  formats retain their previous array mapping and receive no scanner rotation.
- The Cameleon native-array basis is now `[-read,-phase,+normal]`, before
  the patient-coordinate conversion. This agrees with TEST24's empirically
  established export mapping: reverse both in-plane indices, retain slice order.
  Display, contour intersections and DICOM export share this same geometry.
- The exporter derives its direction vectors from the native-array basis and
  the actual stored-pixel flips. It does not apply a second Y/Z conversion.
- 2-D slice-centre spacing is thickness plus gap, including in display and
  contour calculations. The provided SAX example therefore uses 9 mm, not
  70/8 = 8.75 mm. True 3-D data retain FOV/partition spacing.
- Canonical display directions, the orientation compass and registration
  translation labels now use LPS. Imported clinical DICOM contour coordinates
  stay in LPS and do NOT receive the scanner conversion.
- Reconstruction cache signatures carry algorithm version 25. GUI export
  rebuilds stale reconstructions, while direct programmatic export rejects
  stale ImageReconCore objects with instructions to rebuild.
- New exports use a TEST25-specific FrameOfReference UID seed so old and new
  coordinate conventions are not silently labelled as the same frame.
  LastEnhancedDICOMExport diagnostics are excluded from reconstruction keys.

## Verification performed here

- MATLAB syntax parsed successfully for all 10 changed/new MATLAB files using
  MISS_HIT, with UTF-8 input. This is static parsing, not MATLAB execution.
- Independently calculated all eight SAX frame positions from the supplied
  raw prescription and confirmed they match the Y/Z-converted source headers.
- Checked 400 native-voxel/exported-pixel coordinate correspondences, including
  oblique geometry and rectangular matrices, to 1e-10 mm tolerance.
- Independently reconstructed the supplied raw SAX data with unwhitened RSS.
  Mean per-slice correlation against the console, without registration,
  increased from 0.7877 with full offsets to 0.8298 with residual offsets.
  This is a separate numerical check; it is not execution of the MATLAB app
  and does not reproduce its receiver weighting or denoising choices.
- Added MATLAB regression tests for offsets, idempotence, malformed metadata,
  2-D versus 3-D ramps, LPS geometry, slice spacing, oblique voxel mapping and
  an actual exporter/dicomread round trip. Updated existing LPS expectations.

## Remaining validation

MATLAB and CVI were not available here. Run the app and verify anatomical
left/right, anterior/posterior and head/feet directions on the exported images.
The Y/Z conversion is supported by the supplied HFS anatomy and is not a
vendor-certified rule for every scanner/patient-position configuration.
Other patient positions and true 3-D slice ordering need separate validation.

The approximately one-pixel residual discrepancy observed against the console
was NOT corrected with a guessed circshift. The explicit pixel flips and
half-pixel-centred geometry remain paired consistently; any further change to
FFT-origin handling should be established from raw acquisition conventions.

A common patient-coordinate conversion preserves relative geometry between
FCI series. It cannot by itself fix an unrelated within-FCI cross-reference
problem or perform cross-scanner registration. Re-import contours and reset
old manual registration adjustments before assessing the new geometry.

## Use

Follow START_HERE_CVI_TEST.txt. Use a fresh folder and fresh raw-data load.
For a programmatic reconstruction/export:

```matlab
scan = buildimages(scan);  % rebuild from originalcomplexkspace
info = ffc_export_enhanced_mr(scan, 'scan_TEST25.dcm');
```
