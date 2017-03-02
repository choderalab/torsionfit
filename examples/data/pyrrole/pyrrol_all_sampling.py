
from torsionfit import TorsionScanSet as ScanSet
import torsionfit.TorsionFitModel as Model
from torsionfit import sqlite_plus
from pymc import MCMC
from parmed.charmm import CharmmParameterSet
import glob
from pymc import MCMC

structure = 'pyrrol.psf'
stream = 'pyrrol.str'
scan = glob.glob('torsion-scan/*.log')

pyrrol_scan = ScanSet.parse_gauss(scan, structure)
pyrrol_opt = pyrrol_scan.extract_geom_opt()

# set up torsionfit model
param = CharmmParameterSet('../charmm_ff/top_all36_cgenff.rtf', '../charmm_ff/par_all36_cgenff.prm')
model = Model.TorsionFitModelEliminatePhase(param, pyrrol_opt, decouple_n=True, sample_n5=True, stream=stream)

sampler = MCMC(model.pymc_parameters, db=sqlite_plus, dbname='pyrrol_all.db')
sampler.sample(iter=10000)
