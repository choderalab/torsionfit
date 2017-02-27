__author__ = 'Chaya D. Stern'

import simtk.openmm as mm
import simtk.unit as u
import simtk.openmm.app as app

from parmed.charmm import CharmmPsfFile

from pymbar import timeseries

import os
import time
import numpy as np

from torsionfit.utils import logger

# -----------
# Constants
# -----------

kB = u.BOLTZMANN_CONSTANT_kB * u.AVOGADRO_CONSTANT_NA  # Boltzman constant

# ToDo use logger to print out simulation progress
def solvate_system(mol_name, path_to_file):
    pass


def create_openmm_system(mol_name, path_to_files, param):
    """

    :param mol_name:
    :param path_to_files:
    :param param:
    :return:
    """

    # create dict to store simulation data

    psf_vacuum = CharmmPsfFile(os.path.join(path_to_files, "{}.psf".format(mol_name)))
    pdb_vacuum = app.PDBFile(os.path.join(path_to_files, "{}.pdb".format(mol_name)))
    psf_solvate = CharmmPsfFile(os.path.join(path_to_files, "{}_solvated.psf".format(mol_name)))
    pdb_solvate = app.PDBFile(os.path.join(path_to_files, "{}_solvated.pdb".format(mol_name)))

    # create openmm system in solvent
    # Compute box dimensions
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

    psf_solvate.box = (max_crds[0]-min_crds[0],
                 max_crds[1]-min_crds[1],
                 max_crds[2]-min_crds[2], 90.0, 90.0, 90.0)

    system_solvate = psf_solvate.createSystem(param, nonbondedMethod=app.PME,
                                constraints=app.HBonds,
                                nonbondedCutoff = 12.0*u.angstroms,
                                switchDistance=10.0*u.angstroms,)
    positions_solvate = pdb_solvate.getPositions()

    # create openmm system in vacuum

    system_vacuum = psf_vacuum.createSystem(param, nonbondedMethod=app.NoCutoff,
                                            constraints=app.HBonds,
                                            implicitSolvent=None,
                                            removeCMMotion=False)
    positions_vacuum = pdb_vacuum.getPositions()

    # create database to store systems
    database = {'solvated': {'system': system_solvate, 'positions': positions_solvate},
                'vacuum': {'system': system_vacuum, 'positions': positions_vacuum}}
    return database


def generate_simulation_data(database, parameters, solvated=True, n_steps=2500, n_iter=200):
    """

    :param database:
    :param parameters:
    :param solvated:
    :return:
    """

    platform = mm.Platform.getPlatformByName('Reference')

    if solvated:
        system = database['solvated']['system']
        positions = database['solvated']['positions']
    else:
        system = database['vacuum']['system']
        positions = database['vacuum']['positions']

    # create context
    time_step = 2.0*u.femtoseconds
    temperature = 300*u.kelvin
    friction_coef = 1.0/u.picoseconds

    integrator = mm.LangevinIntegrator(temperature, friction_coef, time_step)
    context = mm.Context(system, integrator, platform)

    # set coordinates
    context.setPositions(positions)

    # Minimize
    logger().info('Minimizing Energy...')
    mm.LocalEnergyMinimizer.minimize(context)

    # Simulate, save periodic snapshots
    kT = kB * temperature
    beta = 1.0 / kT

    n_atoms = system.getNumParticles()

    initial_time = time.time()
    x_n = np.zeros([n_iter, n_atoms, 3]) # positions in nm
    u_n = np.zeros([n_iter], np.float64) # energy

    logger.info('Running Simulation...')
    for iteration in range(n_iter):
        logger.info('Iteration {}'.format(iteration))
        integrator.step(n_steps)
        state = context.getState(getEnergy=True, getPositions=True)
        x_n[iteration,:,:] = state.getPositions(asNumpy=True) / u.nanometers
        u_n[iteration] = beta * state.getPotentialEnergy()

    if np.any(np.isnan(u_n)):
        raise Exception("Encountered NaN")

    final_time = time.time()
    elapsed_time = final_time - initial_time
    logger.info('Finished running simulation {} s'.format(elapsed_time))

    # clean up
    del context, integrator

    logger.info('Discarding initial transient equilibration...')
    # Discard initial transient equilibration
    [t0, g, Neff_max] = timeseries.detectEquilibration(u_n)
    x_n = x_n[t0:, :, :]
    u_n = u_n[t0:]

    # subsample to remove correlation
    logger.info('Subsample to remove correlation...')
    indices = time_step.subsampleCorrelatedData(u_n, g=g)
    x_n = x_n[indices, :, :]
    u_n = u_n[indices]

    # store data
    if solvated:
        database['solvated']['x_n'] = x_n
        database['solvated']['u_n'] = u_n
    else:
        database['vacuum']['x_n'] = x_n
        database['vacuum']['u_n'] = u_n

    logger.info("simulation %12.3f s | %5d samples discarded | %5d independent samples remain" % (elapsed_time, t0, len(indices)))

    return database





















