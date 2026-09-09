# SVS relaxometry and 3D multiplanar display — 2 September 2026

## 3D display

- True 3D datasets now expose a **View** selector with axial, coronal and
  sagittal planes.
- Reformatting is performed only after reconstruction. Raw and working
  k-space are not rotated, reflected, interpolated or reordered by the view
  controls.
- The slice selector follows the active anatomical plane.
- Physical read, phase and partition FOVs are carried into each view so
  rectangular voxels retain the correct aspect ratio.
- The orientation compass is calculated from the patient-space basis of the
  displayed plane and follows rotations and reflections. Acquisition-axis
  signs are canonicalised to a radiological display (left=R, right=L; axial
  top=A, bottom=P), and the compass accounts for MATLAB's top-to-bottom image
  row direction.

## Single-voxel FCI

- A scan is recognised as SVS when its sequence or series identifier contains
  `SVS`.
- SVS bypasses spatial FFT reconstruction. For every time/field experiment,
  its fitted signal is `mean(abs(data),1)`, giving an internal magnitude
  array of `[1 1 1 time field]`. The complex mean is retained only for phase
  display and QC. The exact time-by-field arrays are also retained as
  `param.SVS_PROCESSED_COMPLEX_SIGNAL` and
  `param.SVS_PROCESSED_MAGNITUDE_SIGNAL` for inspection.
- **Generate ROI Dispersion Curve** uses that single-voxel signal directly;
  the user is not asked to draw an image ROI.
- Pixelwise relaxation-map generation is disabled for SVS.
- Selecting an SVS acquisition remembers that prescription. Any imaging scan
  subsequently selected in the scan list becomes its localizer and displays
  the remembered voxel. Returning to the SVS uses the most recently selected
  imaging scan; an automatic imaging fallback is used only until one has been
  selected explicitly. The plane nearest the SVS centre is selected
  automatically.
- The voxel outline is the geometric intersection of the prescribed cuboid
  with the displayed localizer plane. Both prescriptions are transformed from
  scanner coordinates into patient coordinates. The explicit
  `OFF_CENTER_FIELD_OF_VIEW_X/Y/Z` planning coordinates are preferred for the
  centre; oriented read/phase/slice offsets remain the compatibility fallback.

## SVS metadata tolerance

- SVS construction no longer dereferences imaging-only fields such as
  `RESOLUTION_FREQUENCY`, `USER_ZERO_FILLING_2D`, `USER_MATRIX_DIMENSION_2D`,
  or `MULTI_PLANAR_EXCITATION` when they are absent.
- SVS preprocessing reshapes the acquired signal directly into
  `[sample phase partition time field receiver]` and bypasses imaging-only
  k-space reorder/padding.

## Verification

### 3 September 2026 constructor hotfix

- Corrected the invocation of defensive Cameleon/header helper methods from
  the `ImageReconCore` constructor. MATLAB parses these helpers as class
  methods, so each call now passes the current object explicitly.
- This resolves the repeated `Undefined function 'getParScalar' for input
  arguments of type 'struct'` failure and applies the same correction to the
  related text, vector, logical, integer and padding helpers.

### 3 September 2026 RS2D spatial-polarity correction

- Corrected the RS2D subject convention to `+X=R`, `+Y=A`, `+Z=H` in the
  anatomical view and compass calculations.
- Separated the scanner prescription basis from the reconstructed-array
  basis. With the Cameleon complex convention used here, reconstructed phase
  and partition indices have the opposite polarity to the prescription
  column and normal stored in `IMAGE_ORIENTATION_SUBJECT`.
- This fixes the measured SVS/localizer mismatch in which the voxel was
  reflected across the localizer centre in the anterior-posterior direction;
  it also selects the correct through-plane partition in a 3D localizer.
- The canonical image pixels remain in the same radiological orientation;
  the correction changes their scanner-coordinate interpretation, the voxel
  overlay and the compass consistently.

### 3 September 2026 simple inversion-recovery fitting

- SVS, imaging ROI dispersion curves and voxelwise relaxation maps now use
  the same unweighted `y(t)=A*exp(-t/T1)+C` model.
- SVS uses `mean(abs(rawdata),1)`; imaging ROIs use the mean magnitude of the
  selected pixels. Magnitude is therefore averaged before fitting rather than
  calculated from a complex average.
- The fitter requires the acquisition to start negative and finish positive,
  but does not hard-code a transition index. Every possible single negative
  prefix is tested and one common index is selected from the aggregate
  scale-normalised SSE across all evolution fields.
- On the supplied 35-field BSA SVS dataset this automatically selects four
  early negative points and reproduces the plausible dispersion trend,
  including the expected quadrupole structure near 50--70 mT.
- Complex phase, phase drift, robust weights, magnitude-noise models,
  multi-field joint fitting and extra physical constraints are no longer in
  the GUI default path. They remain available through
  `T1dispersionAdvanced` and `MapRelaxationAdvanced` if an old analysis needs
  to be reproduced.

Run the MATLAB tests from the project root:

```matlab
results = runtests('tests');
table(results)
```

`test_svs_and_mpr.m` covers magnitude-before-mean SVS reduction, the simple
offset-exponential fit, automatic IR-prefix selection, all three anatomical
plane mappings, lossless display rotation/compass geometry, the supplied
SVS/spin-echo offset regression, and voxel-plane intersection.
`test_reconstruction_geometry.m` continues to guard immutable raw k-space,
centred FFTs, FOV-offset idempotence and receiver whitening.
