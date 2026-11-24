# Enhanced MR DICOM export — TEST19

TEST19 adds a standards-oriented Enhanced MR Image Storage export to the
FFC Processing reconstruction interface.

## Scope of this first implementation

- Exports the complete reconstructed **magnitude** data set as one DICOM Part 10
  Enhanced MR Image Storage SOP Instance.
- Uses Explicit VR Little Endian transfer syntax.
- Preserves the native patient-space reconstruction geometry from
  `ffc_scan_spatial_geometry`; GUI display rotation/mirroring is not baked into
  the DICOM pixels or geometry.
- Stores the five-dimensional FCI magnitude data `[read phase slice time field]`
  using DICOM multi-frame dimensions:
  - **evolution field** -> Stack ID / stack dimension
  - **evolution time** -> Temporal Position Index
  - **slice** -> In-Stack Position Number
- Frame order is `field -> time -> slice`.
- Exact evolution field (mT) and evolution time (ms) are written into each
  frame's Frame Comments so the physical FCI condition is retained even though
  the standard dimension indices are ordinal.
- Pixel data are stored as unsigned 16-bit integers with a global rescale slope
  back to the reconstructed magnitude scale.
- Uses `Content Qualification = RESEARCH`.
- Uses a standards-valid generic whole-body anatomy fallback when a more
  specific anatomical code is not available from the reconstruction metadata.
- No Philips private tags or vendor emulation are used.

## GUI

The **Export Options** panel contains a new **Export Enhanced MR DICOM** button.
Choose a destination `.dcm` file; the scan records basic export metadata under
`scan.param.LastEnhancedDICOMExport` after a successful write.

## Validation status

The existing FFC numerical regression validator is unchanged and should still
pass 22/22 checks. The Enhanced MR writer is statically checked for the required
multi-frame functional-group structure in this build. MATLAB and CVI42 are not
available in the build environment, so the final interoperability test should be
performed by exporting a small reconstructed scan and opening the resulting file
in CVI42 (and/or inspecting it with MATLAB `dicominfo`).
