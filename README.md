# Standard Cartesian RF correction for FFC Processing

The RF subtraction now runs by default in the imaging reconstruction. It supports 2D single-slice and multislice, 3D, multi-receiver, and dynamic Cartesian acquisitions. It uses console acquisition metadata to choose sequential or centric phase order where available and otherwise evaluates candidate orders from the data. Centre-out positive-first and negative-first orders are supported. The scout's sequential back-and-forth transform is evaluated on its stored ky grid. The estimator needs no ADC timestamps.

The RF stage receives the immutable preprocessed complex array, independently examines each acquired receiver, and applies accepted models before the FOV phase ramp, whitening, receiver selection, windowing, FFT and coil combination. For 2D, it validates subtraction independently for each physical slice and frame. For 3D, dimension three remains a kz partition axis, and the previous 3D centric acceptance policy is retained. The three earlier brain fixtures produced byte-identical corrected arrays.

The output keeps the original size and numeric type. For a centre-padded acquisition such as the cine, only measured ky lines are used for estimation and subtraction. The acquired span in that cine was columns 20–109 of the 128-column reconstruction grid. The zero-filled columns stayed zero.

## Installation

1. Back up the active root files `ImageReconCore.m`, `ffc_reconstruct.m` and `ffc_reconstruction_key.m`. The replacements are based on the previously supplied FFC Processing project files.
2. Copy the six production `.m` files from this archive to the active project root. They are `ImageReconCore.m`, `ffc_reconstruct.m`, `ffc_reconstruction_key.m`, `ffc_rf_correction_config.m`, `ffc_correct_rf.m` and `rs2d_predict_coherent_rf.m`. The previous prototype hook and switch are no longer called.
3. Close and reopen the app. Restart an existing MATLAB process pool if active. Run `ffc_setup` from the project root and use `which -all ImageReconCore` and `which -all ffc_reconstruct` to check that the active files are the root copies. The project also contains legacy duplicates.
4. Reopen an acquisition. The updated reconstruction key invalidates older cached images. Warm display and rebuilds using the same source can reuse the cached images or the compact RF model. Forced reconstruction refits from the original k-space. The first RF fit took roughly 3 seconds on the scout and 32–35 seconds on the full multichannel dynamic fixtures in the benchmark environment.

The switch is `enabled = true` in `ffc_rf_correction_config.m`. Set it to `false` for an uncorrected comparison, then reopen the dataset or rebuild. Inspect `scan.rf_correction_info` for the selected order, source phase range, receivers, rejection reasons and fit diagnostics. For 2D, receiver records contain a `slice` array with frame-level outcomes.

## Acceptance and supported data

The model fits outer acquired ky lines and checks predictions in distinct bands closer to the centre. A 2D or sequential fit also requires affirmative corroboration from background margins on central ky lines. This protects against subtracting a tissue feature that resembles the Fourier signature of a stationary RF tone. A fit may be rejected for one receiver or slice while accepted for another. Stationary carriers with abrupt phase changes at the k-space centre remain a limitation, particularly for a near-on-grid peak overlapping anatomy.

The implementation applies to a centred Cartesian grid with at least 32 readout samples and 48 measured phase-encode lines. Explicit non-Cartesian trajectories and spectroscopy or receiver-only acquisitions pass through the RF stage. Acquisition-specific echo trains or variable phase behaviour can also produce an unconfirmed fit.

## Validation

| Series | Data file suffix | Acquired layout | Changed slice-frames | Corrected receivers | Accepted-tone power reduction |
|---|---|---|---:|---:|---:|
| 11604, three-slice FCI | 083757 | 128 × 110 × 3 × 7 × 15 × 13 | 315 of 315 | 2 of 13 | 90.5% |
| 11597, sequential axial scout | 083758 | 128 × 128 × 15 × 1 × 1 × 13 | 10 of 15 | 7 of 13 | 62.0% |
| 11605, cardiac cine | 083741 | 128 × 90 × 9 × 30 × 1 × 13 | 0 of 270 | 0 of 13 | No accepted fits |
| 11606, receiver-only, field off | 083720 | 8192 × 1 × 1 × 1 × 1 × 13 | Receiver-only | 0 of 13 | No phase encoding |

The reduction figures refer to the sum of projected power for **accepted tones** in acquired ky lines excluded from both fitting and acceptance. The corresponding reductions across all discovered candidates were 2.81% for the FCI and 0.12% for the axial scout. Individual visible tones can improve substantially while leaving a small change in total projected power. The real data have no matched clean tissue reference.

The scout's main corrected carrier was reconstructed at 8.409678–8.409682 MHz after applying the header's receiver frequency and effective readout FOV offset. This agrees with the independent field-off spectrum peak at 8.4096804 MHz. The field-off file has no phase encoding and cannot verify central-ky phase prediction. It was used only to check carrier identity.

- [Real images and removed components](validation/real_examples.png)
- [Benchmark summary](validation/summary.csv)
- [Receiver, slice and frame metrics](validation/receiver_slice_frame_metrics.csv)
- [Accepted carrier details](validation/accepted_carriers.csv)

The RF stage and synthetic reconstruction tests ran in GNU Octave 9.4. The synthetic test checks independent slices and receivers, confirmed linear-tone removal, pass-through for an ambiguous compact anatomical feature, zero-pad preservation, whitening and deterministic model reuse. The standard reconstruction integration test uses the FFC Processing functions, with Octave equivalents for MATLAB's `rssq` and the `eye(...,'like',...)` call in the existing whitening function. The real 13-receiver scout also passed the reconstruction-function test with whitening and RSS combination enabled. The MATLAB GUI and process-pool/disk-cache integration were not exercised here.

On the installed app path, run `test_ffc_standard_rf` for the new stage and `test_rs2d_predict_coherent_rf` for the 3D regression. The test and raw fixture reader are optional for routine operation.
