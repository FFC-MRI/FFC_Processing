# FFCProcessing
 Image Processing package for FFC-MRI. Performs preprocessing and reconstruction of FFC-MRI data sets. T1 analysis can be performed using this software or else exported for cluster analysis.

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
<<<<<<< HEAD
  single-channel acquisition. When enabled, all acquired receiver channels are
  prewhitened before the selected channel subset is retained for reconstruction.
- Receiver labels in the multicoil selector match the scanner's zero-based IDs
  (`0` through `N-1`). Internally, those labels map to MATLAB indices `1`
  through `N`.
=======
  single-channel acquisition.
>>>>>>> ab3c741283ada40bce191eb9b97fc9b684810912
- The read, phase and partition dimensions are reconstructed independently;
  rectangular matrices and non-square FOVs are retained rather than forced to
  a square grid.

Run the reconstruction regression tests in MATLAB with:

```matlab
results = runtests('tests/test_reconstruction_geometry.m');
table(results)
```
