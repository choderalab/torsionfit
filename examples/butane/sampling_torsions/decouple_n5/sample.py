import matplotlib.pyplot as plt
import simtk.openmm as mm
from torsionfit import TorsionScanSet as ScanSet
import torsionfit.TorsionFitModel as Model
from torsionfit import sqlite_plus
from pymc import MCMC
from parmed.charmm import CharmmParameterSet
from torsionfit import plots
from matplotlib.backends.backend_pdf import PdfPages

param_to_opt=[('CG331', 'CG321', 'CG321', 'CG331')]
param = CharmmParameterSet('../../data/charmm_ff/top_all36_cgenff.rtf', '../../data/charmm_ff/par_all36_cgenff.prm')
structure = 'structure/butane.psf'
scan = 'MP2_torsion_scan/'
# Print initial guess
print param.dihedral_types[param_to_opt[0]]

# Create a butane scan with torsions on (to compare for difference potential)
#butane_scan = ScanSet.parse_psi4_out(scan, structure)
#butane_scan.compute_energy(param)
#optimized = butane_scan.remove_nonoptimized()
#optimized.compute_energy(param)

# Turn off torsion
param.dihedral_types[('CG331', 'CG321', 'CG321', 'CG331')][1].phi_k=0
param.dihedral_types[('CG331', 'CG321', 'CG321', 'CG331')][0].phi_k=0

# Create butane scan with torsions off
optimized_0 = ScanSet.remove_nonoptimized()
optimized_0.compute_energy(param)

platform = mm.Platform.getPlatformByName('Reference')
model = Model.TorsionFitModelEliminatePhase(param, optimized_0, platform=platform,
                                            param_to_opt=param_to_opt, sample_n5=True, decouple_n=True)
sampler = MCMC(model.pymc_parameters, db=sqlite_plus, dbname='sample.db', verbose=5)

sampler.sample(100000)
