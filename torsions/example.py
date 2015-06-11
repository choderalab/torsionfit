# opemm imports
import simtk.openmm as mm

# ParmEd imports
from parmed.charmm.parameters import CharmmParameterSet

import TorsionScanSet, TorsionFitModel
import glob
from pymc import MCMC

param = CharmmParameterSet('../charmm_ff/top_all36_cgenff.rtf', '../charmm_ff/par_all36_cgenff.prm')
stream = '../structures/Pyrrol/pyrrol.str'
structure = '../structures/Pyrrol/pyrrol.psf'
scan = glob.glob('../structures/Pyrrol/torsion-scan/*.log')
pyrrol_scan = TorsionScanSet.read_scan_logfile(scan, structure)
pyrrol_opt = pyrrol_scan.extract_geom_opt()
#create pymc model
platform = mm.Platform.getPlatformByName('CPU')
model = TorsionFitModel.TorsionFitModel(param, stream, pyrrol_opt, platform=platform)
#update param with missing parameters
model.add_missing(param)
sampler = MCMC(model.pymc_parameters, db='pickle', dbname='pyrrol.database')

sampler.sample(iter=10000)
sampler.db.close()

