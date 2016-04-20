# opemm imports
import simtk.openmm as mm

# ParmEd imports
from parmed.charmm import CharmmParameterSet
#import cProfile
#import pstats
#import StringIO
import torsionfit.TorsionScanSet as TorsionScanSet
import torsionfit.TorsionFitModel as TorsionFitModel
import glob
from pymc import MCMC
import torsionfit.sqlite_plus as db
import time
#from memory_profiler import profile

# Load all parameter, structure and QM log files
param = CharmmParameterSet('data/charmm_ff/top_all36_cgenff.rtf', 'data/charmm_ff/par_all36_cgenff.prm')
streams = ['data/pyrrole/pyrrol.str']#, 'data/pyrrole-2/methyl-pyrrol.str']
structure = ['data/pyrrole/pyrrol.psf']#, 'data/pyrrole-2/methyl-pyrrol.psf']
scans = [glob.glob('data/pyrrole/torsion-scan/*.log')]#, glob.glob('data/pyrrole-2/torsion-scan/*.log')]

# create TorsionScanSet for each fragment
torsion_scans = {}
for mol, scan in zip(structure, scans):
    torsion_scan = TorsionScanSet.read_scan_logfile(scan, mol)
    torsion_scans[mol] = torsion_scan.extract_geom_opt()


#create pymc model
platform = mm.Platform.getPlatformByName('Reference')
model = TorsionFitModel.TorsionFitModel(param, streams, torsion_scans.values(), platform=platform)

sampler = MCMC(model.pymc_parameters, db=db, dbname='test')

#start = time.time()
sampler.sample(iter=10)
#end = time.time()
#print(end-start)

#cProfile.run('sampler.sample(iter=50)', 'statfile')
#sampler.db.close()
#stream = StringIO.StringIO()
#pstats.print_stats()

#print(stream.getvalue())
