# Enhanced MR DICOM export — TEST20

TEST20 extends the Enhanced MR exporter added in TEST19.

## Raw-data demographics

The exporter does not require a console or clinical DICOM. Patient/study metadata
are read from the raw Cameleon dataset:

- top-level `Serie.xml` metadata are retained when the dataset is opened;
- all Cameleon `header.xml` parameters remain available in `scan.param`;
- Patient Name, ID, DOB, sex, study/session date, subject position and acquisition
  timestamp are mapped when present;
- the participant label in `Serie.xml` is used as a PatientID fallback rather than
  writing a generic `FCI` identifier.

`data.dat` itself is the binary complex-sample buffer; identity lives in the raw
XML metadata stored alongside it.

## Participant batch export

The Enhanced DICOM button now offers either the highlighted sequence or all
loaded image sequences for the same participant. Batch export writes one Enhanced
MR object per sequence. All files in the batch deliberately share a
StudyInstanceUID and FrameOfReferenceUID, while each sequence receives its own
SeriesInstanceUID. This avoids the console-writer behaviour where (0020,0052)
FrameOfReferenceUID can vary from series to series for the same participant.

SVS/non-image scans without reconstructed magnitude image data are skipped and
reported in the MATLAB Command Window.
