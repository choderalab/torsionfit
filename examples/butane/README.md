### Butane example

This example goes through all the steps needed to fit the torsion parameters
using Bayesian inference and then propagating the uncertainty to hydration free energy.
Butane is a very simple molecule (only one torsion that does not include any hydrogens) 
so it was used to test this approach and the code. 

1)  Generate the QM torsion scan to fit torsion parameters.
* `torsion_scans` - this folder contains all the relevant scripts and results
* `structure/` - pdb and psf files of butane
* `param/` - edited CHARMM parameter file (non-bonded) used when looking at
energies different parameters contribute.

2)  Sample torsion parameters:
* `c_c_c_c_torsion/` - Only fit the c-c-c-c torsion using RJ, not using RJ, 
omitting n_5 and sampling n_5 (n=5 is usually omitted in CHARMM because it is 
not 'physical'. However, since we end up fitting torsions last, we are essentially
just correcting for everything that was not fit with the other parameters. 
Therefore, the torsion energy is not actually the physical torsion scan we are 
fitting to so it does not have to present something physical
* `all_torsions` - Here I fit all the torsions in butane (5 torsions) because
the c-c-c-c torsion scan should also give us the necessary data for the torsions
involving hydrogen. However, we might have too many parameters vs. data points
(~37 parameters vs 34 data points in the worst case, and 27 parameters in the 
best case)

3)  Propagating uncertainty in hydration free energy.
* `hydration_energy/`