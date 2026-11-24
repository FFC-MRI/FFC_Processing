ISMRM relaxometry colour-map resource

Lipari is used for T1 and R1 maps following:
Fuderer M, Wichtmann B, Crameri F, et al. Color-map recommendation for MR
relaxometry maps. Magn Reson Med. 2025;93(2):490-506. doi:10.1002/mrm.30290.

The recommendation originated with the ISMRM Quantitative MR Study Group and
was approved by the ISMRM Publication Committee and Board. The Lipari palette
is from Fabio Crameri's Scientific Colour Maps and the values here match the
reference resource distributed at https://github.com/mfuderer/colorResources
(associated archive: https://doi.org/10.5281/zenodo.8268884).

FFC Processing applies the recommended logarithm-like transformation to the
colour table, not to the quantitative T1/R1 image values. R1 uses the reversed
Lipari table, matching the MATLAB reference implementation; the paper notes
that consensus was not reached on whether inverse relaxation maps should be
reversed.
