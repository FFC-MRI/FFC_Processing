# cvi42 contour import

This edition can import Circle cvi42 ventricular contours and overlay them on FCI images.

## Workflow

1. Load the FCI dataset normally.
2. Open **CVI Contours > Import workspace...**.
3. Select the `.cvi42ws` (or `.cvi42wsx`) file exported from cvi42.
4. When prompted, select the root folder containing the **clinical DICOM images used to create that cvi42 workspace**. This is required because the workspace contour coordinates are image-relative; DICOM supplies the physical slice orientation, pixel size and position.
5. The registration/display window opens. Choose ED/ES, enable the desired LV/RV boundaries, and adjust X/Y/Z translation and rotation until the contours align with the FCI anatomy.

The importer converts DICOM LPS coordinates to the FCI Processing patient convention (+X right, +Y anterior, +Z head), constructs 3-D surfaces from the serial SAX contours and intersects those surfaces with the currently displayed FCI plane. Therefore it is not limited to identical matrix sizes or identical source/FCI slice numbers.

The initial registration preserves anatomical axes and centres the clinical contour stack on the FCI scan centre. Translation and rotation are then manual because the patient has been repositioned between scanners.

This is a research overlay and must be visually checked before quantitative use.
