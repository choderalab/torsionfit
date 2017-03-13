This directory contains figures of the torsion scans.

### Manifest 
* `psi4_archive/` - Initially, I used the psi4 torsion driver to generate
structures. However, many structures did not converge and even the structures
that did had very high MM energies. Figures and psi4 scripts and output
files of that effort live there.
* `Parameter_energy_contribution.ipynb` - An ipython notebook that looks 
at the mm energy contribution of different parameters
* `Butane_scan_param_energy.pdf` - plot of all parameters contributions
* `Butane_scan_param_energy_mm.pdf` - same as above but does not include
QM energy
* `Butane_total_energy.pdf` - Compares QM and MM torsion scan. 
* `Dihedral_contribution.pdf` - plot of energy contribution for each torsion
in butane. 