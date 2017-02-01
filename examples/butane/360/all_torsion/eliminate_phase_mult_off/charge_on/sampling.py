import simtk.openmm as mm
from torsionfit import TorsionScanSet as ScanSet
import torsionfit.TorsionFitModel as Model
from torsionfit import sqlite_plus
from pymc import MCMC
from parmed.charmm import CharmmParameterSet

param = CharmmParameterSet('../../../../../data/charmm_ff/top_all36_cgenff.rtf',
                           '../../../../../data/charmm_ff/par_all36_cgenff.prm')
structure = '../../../../structure/butane.psf'
scan = '../../../../torsion_scans/DFT_b3lyp/butane_scan_b3lyp_360.log'
butane_scan = ScanSet.parse_psi4(scan, structure)

platform = mm.Platform.getPlatformByName('Reference')

model = Model.TorsionFitModelEliminatePhase(param, butane_scan, platform=platform, decouple_n=True,
                                            param_to_opt=[('CG331', 'CG321', 'CG321', 'CG331'),
                                                          ('HGA3', 'CG331', 'CG321', 'HGA2'),
                                                          ('HGA3', 'CG331', 'CG321', 'CG321'),
                                                          ('HGA2', 'CG321', 'CG321', 'HGA2'),
                                                          ('CG331', 'CG321', 'CG321', 'HGA2')])

sampler = MCMC(model.pymc_parameters, db=sqlite_plus, dbname='butane_360_all_mult_off_charge_on.db', verbose=5)
sampler.sample(100000)
