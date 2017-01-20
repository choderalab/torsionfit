import simtk.openmm as mm
from torsionfit import TorsionScanSet as ScanSet
import torsionfit.TorsionFitModel as Model
from torsionfit import sqlite_plus
from pymc import MCMC
from parmed.charmm import CharmmParameterSet

param = CharmmParameterSet('../../../../../data/charm_ff/top_all36_cgenff.rtf',
                           '../../../../../data/charmm_ff/par_all36_cgenff.prm')
structure = '../../../../structure/butane.psf'
scan = '../../../../torsion_scans/DFT_b3lyp/butane_scan_b3lyp_180.log'
butane_scan = ScanSet.parse_psi4(scan, structure)

platform = mm.Platform.getPlatformByName('Reference')

model = Model.TorsionFitModelEliminatePhase(param, butane_scan, platform=platform, decouple_n=True,
                                            param_to_opt=[('CG331', 'CG321', 'CG321', 'CG331')])

sampler = MCMC(model.pymc_parameters, db=sqlite_plus, dbname='butane_180_mult_off_charge_on.db', verbose=5)
sampler.sample(1000000)
