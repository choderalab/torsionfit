This directory contains figures of the torsion scans.

### Manifest 
* `psi4_archive/` - Initially, I used the psi4 torsion driver to generate
structures. However, many structures did not converge and even the structures
that did had very high MM energies. Figures and psi4 scripts and output
files of that effort live there.
* `distributed_torsion_scan.py` - script to generate torsion scan 
* `B3LYP_torsion_scan/` - directory containing pdb, psi4 input and output files of distributed
qm scan (QM level of theory B3LYP/cc-PVDZ)
* `MP2_torsion_scan/` - directory containing pdb, psi4 input and output files of distributed
qm scan (QM level of theory MP2/6-31G(d)). Psi4 MP2 is 2 orders of magnitude faster than DFT
* `Comparing_QM_methods.ipynb` - ipython notebook that looks at the QM and MM energies of QM minimized structures 
using different QM methods (MP2 vs DFT)
* `butane_b3lyp.png` - b3lyp torsion scan
* `butane_mp2.png` - mp2 torsion scan
* `butane_qm_comparison.png` - comparing energy of both methods
* `Parameter_energy_contribution.ipynb` - An ipython notebook that looks 
at the mm energy contribution of different parameters
* `Butane_scan_param_energy.pdf` - plot of all parameters contributions
* `Butane_scan_param_energy_mm.pdf` - same as above but does not include
QM energy
* `Butane_total_energy.pdf` - Compares QM and MM torsion scan. 
* `Dihedral_contribution.pdf` - plot of energy contribution for each torsion
in butane. 