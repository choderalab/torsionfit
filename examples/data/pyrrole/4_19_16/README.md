### analysis of torsionfit output for pyrrole
This db was generated using pymc sqlite backend (on 4/13/16, analysis was started on 4/19/16) (so state was not saved - can't restart). 
Redundant dihedral bug was fixed but it seems like even less multiplicities got locked in this time. (This happened
because add_missing added missing multiplicity terms before the pymc parameters were initialized so all multiplicities
were turned on. This demonstrates the sensitivity of torsionfit to initial configuration)


###Manifest
* `analysis_4_19_16.ipynb` - ipython notebook where the figures were generated
* `equil_detect_4_19_16.csv` - pymbar.timeseries.detectEquilibrium output for traces of pymc parameters (only torsion parameters k, phase and multiplicity bitstring)
 

