import simtk.openmm as mm
from torsionfit import TorsionScanSet as ScanSet
import torsionfit.TorsionFitModel as Model
from torsionfit import sqlite_plus
from pymc import MCMC
from parmed.charmm import CharmmParameterSet

param = CharmmParameterSet('../data/charmm_ff/top_all36_cgenff.rtf', '../data/charmm_ff/par_all36_cgenff.prm')
structure = 'butane.psf'
scan = 'butane_scan_b3lyp_4.log'
butane_scan = ScanSet.parse_psi4(scan, structure)

platform = mm.Platform.getPlatformByName('Reference')

model = Model.TorsionFitModelEliminatePhase(param, butane_scan, platform=platform,
                                            param_to_opt=[('CG331', 'CG321', 'CG321', 'CG331'),
                                                          ('HGA3', 'CG331', 'CG321', 'HGA2'),
                                                          ('HGA3', 'CG331', 'CG321', 'CG321'),
                                                          ('HGA2', 'CG321', 'CG321', 'HGA2'),
                                                          ('CG331', 'CG321', 'CG321', 'HGA2')])

sampler = MCMC(model.pymc_parameters, db=sqlite_plus, dbname='butane_all_torsions', verbose=5)
sampler.sample(100000)
