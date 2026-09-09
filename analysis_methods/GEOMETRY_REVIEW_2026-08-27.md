# FOV offset and rotation review — 27 August 2026

## Reconstruction order

Every image build now follows one deterministic sequence:

1. Copy immutable, preprocessed acquisition data from
   `originalcomplexkspace` into the working `complexkspace`.
2. Reset the working-copy offset flag.
3. Apply the scanner's 1D/2D FOV-offset phase ramps once. Apply the 3D ramp
   only for a true 3D acquisition; the third dimension of a 2D multislice scan
   is a slice label and receives no kz ramp.
4. Prewhiten the complete acquired receiver array when selected, then retain
   the chosen receiver subset.
5. Window, zero-fill and Fourier transform the spatial dimensions.
6. Combine receivers and generate magnitude/phase images.
7. Apply the accumulated GUI rotation/reflection to reconstructed arrays using
   exact dimension permutation and reversal operations.

## Cause of the rebuild shift

The previous rotate and flip callbacks used `imwarp` on
`originalcomplexkspace`. This changed the source used by later builds. The
header offsets still described the original read, phase and partition axes, so
the next build applied their phase ramps to a differently oriented and
interpolated k-space array. With non-zero offsets this appeared as an image
translation after rebuilding; rectangular matrices made the inconsistency more
obvious.

## Current guarantees

- `originalcomplexkspace` is never modified by display rotation or reflection.
- FOV-offset application is explicitly idempotent for each working copy.
- A rebuild always starts from the same unshifted preprocessed data.
- The same accumulated `geomT` controls reconstructed-array orientation, the
  orientation compass and whether read/phase display extents are swapped.
- Display `XData`/`YData` now describe pixel centres and explicit
  `XLim`/`YLim` describe the full physical FOV edges. This prevents MATLAB's
  automatic half-pixel expansion from changing the aspect ratio when a
  rectangular matrix is rotated.
- Ninety-degree rotations and flips use no interpolation and preserve numeric
  values exactly.
- `fov_offset_info` records the offsets, FOVs and 2D/3D decision used by the
  most recent build.
