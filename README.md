# FFC Processing — Performance Edition

From MATLAB, make this extracted folder current and run:

```matlab
ffc_start
```

The active app is `FFCProcessing.m`. Start with [PERFORMANCE_EDITION.md](PERFORMANCE_EDITION.md) for changes, validation status, tests and benchmarking. The original App Designer file is retained in `legacy_app/FFCProcessing_original.mlapp`.

This edition includes reconstruction reuse, bounded FFT/whitening intermediates, persistent image displays, cached filter previews and managed background work. MATLAB runtime tests and real-data timings have not been run in the preparation environment.

 Image Processing package for FFC-MRI. Performs preprocessing and reconstruction of FFC-MRI data sets. T1 analysis can be performed using this software or else exported for cluster analysis.

## SVS and 3D display

- Axial, coronal and sagittal display of true 3D reconstructions, with a
  radiological patient-orientation-aware compass and physical FOV scaling.
  Acquisition-axis signs are canonicalised so axial views consistently show
  R on the left and A at the top, including when 2D and 3D headers use
  opposite row/column signs.
- Localised FCI SVS support. For each experiment the fitted voxel signal is
  `mean(abs(data),1)`, matching direct manual extraction. The complex mean is
  retained only for phase display/QC. SVS can be analysed with **Generate ROI
  Dispersion Curve** without drawing an ROI.
- After selecting an SVS acquisition, selecting any imaging dataset in the
  scan list displays that image with the selected SVS prescription overlaid.

## Inversion-recovery fitting

- SVS, imaging ROIs and relaxation maps use the same ordinary unweighted
  model: `y(t)=A*exp(-t/T1)+C`.
- Magnitude is averaged first: over sampled points for SVS and over selected
  pixels for an imaging ROI.
- The acquisition is treated as inversion recovery, with one early-negative
  to late-positive transition. Its location is selected from the data across
  all fields; no fixed number of negative points is hard-coded.
- The default fit does not use complex phase, phase-drift correction, robust
  weights, Rician/noise-floor modelling, joint-field fitting or additional
  physical constraints. The former implementations remain callable as
  `T1dispersionAdvanced` and `MapRelaxationAdvanced` for reproducibility.

See `SVS_MPR_UPDATE_2026-09-02.md` for behaviour and validation details.

Needs the following MATLAB toolboxes:

Image Processing Toolbox  
Optimization Toolbox  
Curve Fitting Toolbox  
Signal Processing Toolbox  
Parallel Computing Toolbox  

## Reconstruction notes

- Receive-noise whitening is controlled from the **Noise whitening** checkbox
  in **Reconstruction Options**. It is available for multichannel data,
  defaults on when more than one receiver is present, and is forced off for a
  single-channel acquisition. When enabled, all acquired receiver channels are
  prewhitened before the selected channel subset is retained for reconstruction.
- Receiver labels in the multicoil selector match the scanner's zero-based IDs
  (`0` through `N-1`). Internally, those labels map to MATLAB indices `1`
  through `N`.
- The read, phase and partition dimensions are reconstructed independently;
  rectangular matrices and non-square FOVs are retained rather than forced to
  a square grid.
- Prescribed FOV offsets are applied once to a fresh working copy of the
  immutable preprocessed k-space on every rebuild. GUI rotations and flips are
  exact image-plane operations applied after reconstruction; they never warp
  or overwrite acquisition k-space.
- Image axes use explicit physical FOV edges and pixel-centre coordinates, so
  rotating a rectangular matrix swaps its read/phase extents without changing
  its physical aspect ratio.

Run all regression tests in MATLAB with:

```matlab
results = runtests('tests');
table(results)
```
