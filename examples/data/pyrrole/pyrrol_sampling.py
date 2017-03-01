
from torsionfit import TorsionScanSet as ScanSet
import torsionfit.TorsionFitModel as Model
from torsionfit import sqlite_plus
from pymc import MCMC
from parmed.charmm import CharmmParameterSet
import glob
from pymc import MCMC

structure = 'pyrrol.psf'
scan = glob.glob('torsion-scan/*.log')

pyrrol_scan = ScanSet.parse_gauss(scan, structure)
pyrrol_opt = pyrrol_scan.extract_geom_opt()

# set up torsionfit model
param_to_optimize = [('CG331', 'CG321', 'CG321', 'NG3C51'),
                    ('CG321', 'CG321', 'NG3C51', 'CG2R51'),
                    ('CG321', 'CG321', 'NG3C51', 'CG251O'),
                    ('CG2D2', 'CG251O', 'NG3C51', 'CG321'),]
param = CharmmParameterSet('../charmm_ff/top_all36_cgenff.rtf', '../charmm_ff/par_all36_cgenff.prm', 'pyrrol.str')
model = Model.TorsionFitModelEliminatePhase(param, pyrrol_opt, decouple_n=True, sample_n5=True,
                                            param_to_opt=param_to_optimize )

sampler = MCMC(model.pymc_parameters, db=sqlite_plus, dbname='pyrrol_4.db')
sampler.sample(iter=10000)