# Enhanced DICOM TEST 22

- Corrects RS2D subject-coordinate RAH (+R,+A,+H) to DICOM patient-coordinate LPS (+L,+P,+H).
- Stores frames with the corresponding 180-degree in-plane reordering so anatomy agrees with conventional DICOM orientation labels.
- Re-reads raw XML identity at export time and searches participant/study XML context; no console DICOM dependency.
- Removes invented `FCI` Patient ID fallback: missing Type-2 identity fields are left empty rather than wrong.
- Uses raw Study ID/description and raw Series ID where available.
- Separates Study Date/Time from Series/Acquisition/Content Date/Time.
- Writes Patient Weight when present in raw metadata.
