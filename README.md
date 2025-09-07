# MinAlignConf

Requirements:
(a) Rdkit
(b) Numpy
(c)Pandas
(d) tKinter

This GUI is designed for sdf file generation from .csv file with SMILES notations. The structures are geometrically optimized and aligned. It also generates conformations.

The starting .csv file should contain three columns: (a) First column contains SMILES notations of structures with header exactly as 'SMILES'; (b) Second column contains name of compounds with a header (no specification); (c) third columns contains biological activity (i.e., pIC50, pKi) with a header (no specification).

The geometrical optimization of structures is done with UFF (Universal force field). The minimized structures are automatically downloaded as '*_minimized.sdf'.

The alignment and conformer generation are done with Rdkit. The aligned conformations and multiple conformations are also saved automatically.


##Citation:
Banerjee, T., Mitra, S., Sar, S., Halder, A. K., Panda, P., Ghosh, N. Structural insights into molecular and cellular level FXR binding potentials of GW4064 and LY2562175 hybrids by multi in silico modelling analyses. Journal of Molecular Modeling. 2025, 31, 119. https://doi.org/10.1007/s00894-025-06336-5
