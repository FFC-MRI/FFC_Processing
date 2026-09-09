# CVI contour overlay test build 2 — 2026-09-08

This test build adds import and FCI-plane projection of Circle cvi42 short-axis ventricular contours.

## Added

- `CVI Contours` menu in `FFCProcessing`.
- Binary `.cvi42ws` and XML `.cvi42wsx` readers.
- Matching of cvi42 image/frame keys to the original clinical DICOM series.
- Conversion from source-image contour pixels to DICOM patient coordinates and then to the FCI patient-coordinate convention.
- ED/ES labelling when exactly two LV endocardial contour phases are present.
- 3-D open surfaces for LV endocardium, LV epicardium and RV endocardium.
- Intersection of those surfaces with the currently displayed FCI plane.
- Manual 6-degree-of-freedom rigid registration controls.
- `cvi_test_import` command-line diagnostic.

## Import in the app

1. Start the suite with `ffc_start` and load/reconstruct the FCI scan normally.
2. Choose **CVI Contours > Import workspace...**.
3. Select the cvi42 `.cvi42ws` or `.cvi42wsx` file.
4. Select the root folder containing the original clinical DICOM study used by cvi42.
5. The registration window opens. The source group with the most matched ventricular contours is selected automatically, and ED is preferred when identifiable.
6. Adjust X/Y/Z translation and X/Y/Z rotation until the projected contours coincide with the FCI anatomy.

The zero registration deliberately centres the clinical contour stack on the FCI prescription centre. It does **not** assume that scanner coordinates are shared across scanners.

## Test-study expectation

The supplied test workspace was independently decoded during development as 61 supported ventricular contour objects: 21 LV endocardial, 22 LV epicardial and 18 RV endocardial. The principal short-axis cine group contains 55 of these contours and should be selected automatically as the default source group.

## Limitations of this first build

- Registration is manual; no landmark or automatic multimodal registration yet.
- Only SAX LV endocardial, LV epicardial and RV endocardial cvi42 contour names are imported.
- The clinical DICOMs are needed at import time because the cvi42 workspace points alone do not establish the complete source slice position/orientation in patient space.
- The overlay is intended for research processing and should be visually checked before quantitative use.
- This package has been statically reviewed and its supplied workspace format independently decoded, but MATLAB itself was not available in the build environment. If MATLAB reports an error, preserve the full error/stack and the test build can be corrected directly.


## Test build 2 launcher / UI fix

The main window now has two explicit buttons at the bottom of the Scan List panel:
`Import CVI contours` and `CVI registration`.  Use `ffc_start_cvi` to launch this build.
The launcher resolves `FFCProcessingCVI.m` from this folder and prints the path, which avoids accidentally opening an older cached/path copy or the legacy `.mlapp`.


## TEST 4 launcher fix
The active class is uniquely named `FFCProcessingCVI`, so a stale `FFCProcessing.mlapp` left behind when extracting over an older installation cannot take precedence in MATLAB class resolution.


## TEST 6 decompression compatibility fix

TEST 4 used Java InflaterInputStream.read(buffer,...) with a MATLAB numeric
array when zlibdecode was unavailable. On some MATLAB/JVM combinations Java
writes into a temporary copy, so MATLAB received a zero-filled output buffer.
TEST 6 keeps bulk copying within Java when possible and has a scalar Java
read fallback that is compatible with older MATLAB releases.


TEST 9: Added live R/L, A/P and H/F translation sliders (numeric fields retained).
