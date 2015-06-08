# opemm imports
import simtk.unit as u
import simtk.openmm as mm
import simtk.openmm.app as app

# ParmEd imports
from chemistry.charmm import  CharmmPsfFile
from chemistry.charmm.parameters import CharmmParameterSet

import TorsionScanSet, TorsionFitModel
import glob
from pymc import MCMC
import matplotlib
import matplotlib.pyplot as plt
import time
import cProfile
import pstats
import StringIO

param = CharmmParameterSet('../charmm_ff/top_all36_cgenff.rtf', '../charmm_ff/par_all36_cgenff.prm')
stream = '../structures/Pyrrol/pyrrol.str'
structure = '../structures/Pyrrol/pyrrol.psf'
scan = glob.glob('../structures/Pyrrol/torsion-scan/*.log')
pyrrol_scan = TorsionScanSet.read_scan_logfile(scan, structure)
pyrrol_opt = pyrrol_scan.extract_geom_opt()
#create pymc model
platform = mm.Platform.getPlatformByName('Reference')
model = TorsionFitModel.TorsionFitModel(param, stream, pyrrol_opt, platform=platform)
#update param with missing parameters
model.add_missing(param)
sampler = MCMC(model.pymc_parameters)

cProfile.run('sampler.sample(iter=10, burn=0, thin=1)', 'statfile')
stream = StringIO.StringIO()
stats = pstats.Stats('statfile', stream=stream)
stats.print_stats()

print(stream.getvalue())
