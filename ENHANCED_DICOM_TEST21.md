# Enhanced DICOM TEST 21

- Preserves the folder originally selected when raw data are loaded recursively (`openRootPath`).
- Defaults DICOM exports to `<original open folder>/dicom`.
- Creates that `dicom` folder automatically when possible.
- Selected-series `uiputfile` opens in that folder.
- Participant batch `uigetdir` opens in that folder.
- Falls back to the individual dataset folder for scans loaded before TEST 21, then to MATLAB `pwd`.
