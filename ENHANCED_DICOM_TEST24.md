# Enhanced DICOM TEST 24

TEST24 corrects the image-to-patient mapping using an empirical comparison of the
supplied raw `data.dat` (3-slice FCI acquisition) with the matching classic
console DICOM series.

## Measured mapping

- `IMAGE_ORIENTATION_SUBJECT` is retained directly as DICOM
  `ImageOrientationPatient`.
- `OFF_CENTER_FIELD_OF_VIEW_X/Y/Z` is retained directly as the stack centre in
  patient coordinates.
- Reconstructed MATLAB slice index maps directly to the corresponding DICOM
  slice index; **slice order is not reversed**.
- Before MATLAB column-major pixel data are written as DICOM row-major data,
  **both in-plane MATLAB dimensions are reversed**. This reproduces the console
  pixel orientation and keeps anatomy consistent with the DICOM orientation
  labels.
- Per-slice `ImagePositionPatient` is generated from the stack centre, direct
  prescription normal, and the 2-D centre spacing (`slice thickness + gap`).

The change is confined to Enhanced-DICOM export and does not alter reconstruction
or CVI contour registration inside MATLAB.
