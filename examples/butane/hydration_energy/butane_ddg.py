from torsionfit.hydration_energy import energy
from parmed.charmm import CharmmParameterSet
from parmed.topologyobjects import DihedralType
import cPickle as pickle
from torsionfit import sqlite_plus, utils
import simtk.openmm as mm
import simtk.unit as u
import simtk.openmm.app as app
from pymbar import timeseries
import numpy as np
from tqdm import *

param = CharmmParameterSet('../../data/charmm_ff/par_all36_cgenff.prm', '../../data/charmm_ff/top_all36_prot.rtf',
                          '../../data/charmm_ff/toppar_water_ions.str')
database = energy.create_openmm_system('butane', '../structure/', param)

solvated = pickle.load(open('solvated.pickle', 'r'))
vacuum = pickle.load(open('vacuum_db.pickle', 'r'))

database['solvated']['x_n'] = solvated['solvated']['x_n']
database['solvated']['u_n'] = solvated['solvated']['u_n']
database['vacuum']['x_n'] = vacuum['vacuum']['x_n']
database['vacuum']['u_n'] = vacuum['vacuum']['u_n']

# load sampled parameters
db = sqlite_plus.load('../butane_n5_decouple_n.db')

# get parameters
param_to_opt = utils.get_sampled_torsions(db)
equil_t = np.zeros([6,3])
for i in range(6):
    torsion_k = '{}_{}_K'.format(param_to_opt[0], str(i+1))
    equil_t[i] = timeseries.detectEquilibration(db.trace(torsion_k)[:])

min_n = max(equil_t[:,0])
# draw indices randomly to
idx = np.random.randint(min_n, len(db.trace(torsion_k)[:]), 200)

# Add missing parameters
multiplicities = [1, 2, 3, 4, 5, 6]
p = ('CG331', 'CG321', 'CG321', 'CG331')
per = []
for i in range(len(param.dihedral_types[p])):
    per.append(param.dihedral_types[p][i].per)
for j in multiplicities:
    if j not in per:
        param.dihedral_types[p].append(DihedralType(0, j, 0))


# parameterize without deep copy parameters
psf_solvate = database['solvated']['structure']
psf_solvate.load_parameters(param, copy_parameters=False)
system_solvate = psf_solvate.createSystem(nonbondedMethod=app.PME,
                                          constraints=app.HBonds,
                                          nonbondedCutoff = 12.0*u.angstroms,
                                          switchDistance=10.0*u.angstroms)
database['solvated']['system'] = system_solvate

# create openmm system in vacuum
psf_vacuum = database['vacuum']['structure']
psf_vacuum.load_parameters(param, copy_parameters=False)
system_vacuum = psf_vacuum.createSystem(nonbondedMethod=app.NoCutoff,
                                            constraints=app.HBonds,
                                            implicitSolvent=None,
                                            removeCMMotion=False)
database['vacuum']['system'] = system_vacuum
# create context and add to database (this will be done inside energy)
 # create context
time_step = 2.0*u.femtoseconds
temperature = 300*u.kelvin
friction_coef = 1.0/u.picoseconds
integrator = mm.LangevinIntegrator(temperature, friction_coef, time_step)
integrator_2 = mm.LangevinIntegrator(temperature, friction_coef, time_step)

platform = mm.Platform.getPlatformByName('Reference')

context_solvated = mm.Context(database['solvated']['system'], integrator, platform)
database['solvated']['context'] = context_solvated
context_vacuum = mm.Context(database['vacuum']['system'], integrator_2, platform)
database['vacuum']['context'] = context_vacuum

# Test that I get zero when parameters are equal
energy.zwanzig(database['vacuum'])

# Set phase to 0
param.dihedral_types[p][1].phase = 0.0

ddg_vacuum = np.zeros(len(idx))
ddg_solvated = np.zeros(len(idx))
for i in tqdm(range(len(idx))):
    utils.param_from_db(param, db, idx[i], decouple_n=True)
    u_vac = energy.zwanzig(database['vacuum'])
    ddg_vacuum[i] = u_vac / u_vac.unit
    u_solv = energy.zwanzig(database['solvated'])
    ddg_solvated[i] = u_solv / u_solv.unit

pickle.dump(ddg_vacuum, open('ddg_vacum.pickle', 'w'))
pickle.dump(ddg_solvated, open('ddg_solv.pickle', 'w'))
