# TEST 17 — ISMRM relaxometry colour maps

This build replaces the generic MATLAB colour table used for T1/R1 maps with the
ISMRM Quantitative MR Study Group recommendation from Fuderer et al. (MRM 2025).

- T1: Lipari, logarithm-like colour-table remapping (default).
- R1: reversed Lipari, logarithm-like colour-table remapping, matching the
  authors' MATLAB reference implementation.
- The T1/R1 numeric image is **not** log-transformed. The logarithm-like mapping
  is applied to the colour table, so the colour bar remains in physical units.
- Zero/non-finite map values display as black (invalid).
- A labelled colour bar is shown for relaxometry maps.
- TIFF relaxometry exports use the same Lipari mapping as the GUI.
- The existing map-log state button is repurposed: default `ISMRM log-Lipari`;
  switching it off gives `Linear Lipari` for diagnostic comparison.
- Includes the older-MATLAB `uitable` string-to-char compatibility fix from the
  TEST 16 results-table hotfix.

Reference: Fuderer M, Wichtmann B, Crameri F, et al. Color-map recommendation
for MR relaxometry maps. Magn Reson Med. 2025;93(2):490-506.
doi:10.1002/mrm.30290. Palette/reference resources:
https://doi.org/10.5281/zenodo.8268884
