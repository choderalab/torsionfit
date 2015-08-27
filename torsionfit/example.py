# opemm imports
import simtk.openmm as mm

# ParmEd imports
from parmed.charmm import CharmmParameterSet
#import cProfile
#import pstats
#import StringIO
import TorsionScanSet, TorsionFitModel
import glob
from pymc import MCMC
#from memory_profiler import profile

param = CharmmParameterSet('../charmm_ff/top_all36_cgenff.rtf', '../charmm_ff/par_all36_cgenff.prm')
stream = '../examples/pyrrole/pyrrol.str'
structure = '../examples/pyrrole/pyrrol.psf'
scan = glob.glob('../examples/pyrrole/torsion-scan/*.log')
pyrrol_scan = TorsionScanSet.read_scan_logfile(scan, structure)
pyrrol_opt = pyrrol_scan.extract_geom_opt()
#create pymc model
platform = mm.Platform.getPlatformByName('Reference')
model = TorsionFitModel.TorsionFitModel(param, stream, pyrrol_opt, platform=platform)
#update param with missing parameters
model.add_missing(param)
sampler = MCMC(model.pymc_parameters, db='sqlite', dbname='pyrrol.database', verbose=5)

sampler.sample(iter=10000)

#cProfile.run('sampler.sample(iter=50)', 'statfile')
#sampler.db.close()
#stream = StringIO.StringIO()
#pstats.print_stats()

#print(stream.getvalue())
