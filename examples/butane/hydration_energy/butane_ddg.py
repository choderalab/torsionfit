from torsionfit.hydration_energy import energy
from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed.topologyobjects import DihedralType
import cPickle as pickle
from torsionfit import sqlite_plus, utils
import simtk.openmm as mm
import simtk.unit as u
import simtk.openmm.app as app
from pymbar import timeseries
import numpy as np
from tqdm import *

param_1 = CharmmParameterSet('../../data/charmm_ff/par_all36_cgenff.prm', '../../data/charmm_ff/top_all36_prot.rtf',
                          '../../data/charmm_ff/toppar_water_ions.str')
param_2 = CharmmParameterSet('../../data/charmm_ff/par_all36_cgenff.prm', '../../data/charmm_ff/top_all36_prot.rtf',
                          '../../data/charmm_ff/toppar_water_ions.str')
database = energy.create_openmm_system('butane', '../structure/', param_1)

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
idx_1 = np.random.randint(min_n, len(db.trace(torsion_k)[:]), 200)
idx_2 = np.random.randint(min_n, len(db.trace(torsion_k)[:]), 200)
#samples_with_replacement = [[data[i] for i in index] for index in idx]

# Add missing parameters
multiplicities = [1, 2, 3, 4, 5, 6]
p = ('CG331', 'CG321', 'CG321', 'CG331')
per = []
for i in range(len(param_1.dihedral_types[p])):
    per.append(param_1.dihedral_types[p][i].per)
for j in multiplicities:
    if j not in per:
        param_1.dihedral_types[p].append(DihedralType(0, j, 0))
        param_2.dihedral_types[p].append(DihedralType(0, j, 0))

# Set phase to 0
param_1.dihedral_types[p][1].phase = 0.0
param_2.dihedral_types[p][1].phase = 0.0


# Create 2 systems tied to the 2 different parameter sets for solvated and vacuum
# parameterize without deep copy parameters

pdb_solvate = app.PDBFile('../structure/butane_solvated.pdb')

coords = pdb_solvate.positions
min_crds = [coords[0][0], coords[0][1], coords[0][2]]
max_crds = [coords[0][0], coords[0][1], coords[0][2]]

for coord in coords:
    min_crds[0] = min(min_crds[0], coord[0])
    min_crds[1] = min(min_crds[1], coord[1])
    min_crds[2] = min(min_crds[2], coord[2])
    max_crds[0] = max(max_crds[0], coord[0])
    max_crds[1] = max(max_crds[1], coord[1])
    max_crds[2] = max(max_crds[2], coord[2])

psf_solvate_1 = CharmmPsfFile('../structure/butane_solvated.psf')
psf_solvate_1.box = (max_crds[0]-min_crds[0],
                     max_crds[1]-min_crds[1],
                     max_crds[2]-min_crds[2], 90.0, 90.0, 90.0)
database['solvated']['structure_1'] = psf_solvate_1
psf_solvate_1.load_parameters(param_1, copy_parameters=False)
system_solvate_1 = psf_solvate_1.createSystem(nonbondedMethod=app.PME,
                                          constraints=app.HBonds,
                                          nonbondedCutoff = 12.0*u.angstroms,
                                          switchDistance=10.0*u.angstroms)
database['solvated']['system_1'] = system_solvate_1

psf_solvate_2 = CharmmPsfFile('../structure/butane_solvated.psf')
psf_solvate_2.box = (max_crds[0]-min_crds[0],
                     max_crds[1]-min_crds[1],
                     max_crds[2]-min_crds[2], 90.0, 90.0, 90.0)
database['solvated']['structure_2'] = psf_solvate_2
psf_solvate_2.load_parameters(param_2, copy_parameters=False)
system_solvate_2 = psf_solvate_2.createSystem(nonbondedMethod=app.PME,
                                          constraints=app.HBonds,
                                          nonbondedCutoff = 12.0*u.angstroms,
                                          switchDistance=10.0*u.angstroms)
database['solvated']['system_2'] = system_solvate_2

# create openmm system in vacuum
psf_vacuum_1 = CharmmPsfFile('../structure/butane.psf')
database['vacuum']['structure_1'] = psf_vacuum_1
psf_vacuum_1.load_parameters(param_1, copy_parameters=False)
system_vacuum_1 = psf_vacuum_1.createSystem(nonbondedMethod=app.NoCutoff,
                                            constraints=app.HBonds,
                                            implicitSolvent=None,
                                            removeCMMotion=False)
database['vacuum']['system_1'] = system_vacuum_1

psf_vacuum_2 = CharmmPsfFile('../structure/butane.psf')
database['vacuum']['structure_2'] = psf_vacuum_2
psf_vacuum_2.load_parameters(param_2, copy_parameters=False)
system_vacuum_2 = psf_vacuum_2.createSystem(nonbondedMethod=app.NoCutoff,
                                            constraints=app.HBonds,
                                            implicitSolvent=None,
                                            removeCMMotion=False)
database['vacuum']['system_2'] = system_vacuum_2
# create context and add to database (this will be done inside energy)
 # create context
time_step = 2.0*u.femtoseconds
temperature = 300*u.kelvin
friction_coef = 1.0/u.picoseconds
integrator_1 = mm.LangevinIntegrator(temperature, friction_coef, time_step)
integrator_2 = mm.LangevinIntegrator(temperature, friction_coef, time_step)
integrator_3 = mm.LangevinIntegrator(temperature, friction_coef, time_step)
integrator_4 = mm.LangevinIntegrator(temperature, friction_coef, time_step)

platform = mm.Platform.getPlatformByName('CUDA')

context_solvated_1 = mm.Context(database['solvated']['system_1'], integrator_1, platform)
database['solvated']['context_1'] = context_solvated_1
context_solvated_2 = mm.Context(database['solvated']['system_2'], integrator_2, platform)
database['solvated']['context_2'] = context_solvated_2
context_vacuum_1 = mm.Context(database['vacuum']['system_1'], integrator_3, platform)
database['vacuum']['context_1'] = context_vacuum_1
context_vacuum_2 = mm.Context(database['vacuum']['system_2'], integrator_4, platform)
database['vacuum']['context_2'] = context_vacuum_2



dg_vacuum = np.zeros(len(idx_1))
dg_solvated = np.zeros(len(idx_1))
for i in tqdm(range(len(idx_1))):
    utils.param_from_db(param_1, db, idx_1[i], decouple_n=True)
    utils.param_from_db(param_2, db, idx_2[i], decouple_n=True)
    u_vac = energy.zwanzig(database['vacuum'])
    dg_vacuum[i] = u_vac / u_vac.unit
    u_solv = energy.zwanzig(database['solvated'])
    dg_solvated[i] = u_solv / u_solv.unit

pickle.dump(dg_vacuum, open('dg_vacum.pickle', 'w'))
pickle.dump(dg_solvated, open('dg_solv.pickle', 'w'))
