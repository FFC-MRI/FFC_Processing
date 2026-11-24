# Enhanced DICOM TEST 23

## Geometry correction

- Fixes oblique Enhanced-MR slice orientation.
- `IMAGE_ORIENTATION_SUBJECT` is exported directly as DICOM `ImageOrientationPatient`; no RAH-to-LPS sign reflection is applied.
- This is based on the supplied console DICOMs, where the private/raw `IMAGE_ORIENTATION_SUBJECT` value exactly matches `(0020,0037)` for axial, 2CH, SAX, 4CH, true-SAX, FCI and cine series.
- Native FCI reconstruction reverses prescription phase and slice/partition index directions. Export therefore flips the native phase dimension and reverses source slice order while leaving read order unchanged.
- DICOM slice positions advance along `cross(IOP_row, IOP_column)`, preserving the correct 3-D relationship between SAX and long-axis prescriptions.
- Patient/study metadata and batch export behaviour from TEST 22 are retained.
