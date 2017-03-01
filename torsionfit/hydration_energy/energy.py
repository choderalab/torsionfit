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
from torsionfit.TorsionScanSet import TorsionScanSet as ScanSet

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
    # ToDo make sure to add missing parameters before creating the system

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

    # parameterize without deep copy parameters
    psf_solvate.load_parameters(param, copy_parameters=False)
    system_solvate = psf_solvate.createSystem(nonbondedMethod=app.PME,
                                              constraints=app.HBonds,
                                              nonbondedCutoff = 12.0*u.angstroms,
                                              switchDistance=10.0*u.angstroms)
    positions_solvate = pdb_solvate.getPositions()

    # create openmm system in vacuum
    psf_vacuum.load_parameters(param, copy_parameters=False)
    system_vacuum = psf_vacuum.createSystem(nonbondedMethod=app.NoCutoff,
                                            constraints=app.HBonds,
                                            implicitSolvent=None,
                                            removeCMMotion=False)
    positions_vacuum = pdb_vacuum.getPositions()

    # create database to store systems
    database = {'solvated': {'structure': psf_solvate, 'system': system_solvate, 'positions': positions_solvate},
                'vacuum': {'structure': psf_vacuum, 'system': system_vacuum, 'positions': positions_vacuum}}
    #return hydration_endpoints(database)
    return database

# class hydration_endpoints():
#
#     def __init__(self, database):
#         self.vacuum = ScanSet(positions=database['vacuum']['positions'], topology=)
#         self.solvated = ScanSet()


def generate_simulation_data(database, parameters, solvated=True, n_steps=2500, n_iter=200):
    """

    :param database:
    :param parameters:
    :param solvated:
    :return:
    """

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
    platform = _determine_fastest_platform()
    platform_name = platform.getName()
    logger().info('Using {} platform'.format(platform_name))
    if platform_name == 'CUDA':
        prop = dict(CudaPrecision='mixed')
    elif platform_name == 'OpenCL':
        prop = dict(OpenCLPrecision='mixed')
    elif platform_name == 'CPU':
        platform = mm.Platform.getPlatformByName('Reference')
        prop = None

    if prop is None:
        context = mm.Context(system, integrator, platform)
    else:
        context = mm.Context(system, integrator, platform, prop)

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

    logger().info('Running Simulation...')
    for iteration in range(n_iter):
        logger().info('Iteration {} in {} s'.format(iteration, time.time() - initial_time))
        integrator.step(n_steps)
        state = context.getState(getEnergy=True, getPositions=True)
        x_n[iteration,:,:] = state.getPositions(asNumpy=True) / u.nanometers
        u_n[iteration] = beta * state.getPotentialEnergy()

    if np.any(np.isnan(u_n)):
        raise Exception("Encountered NaN")

    final_time = time.time()
    elapsed_time = final_time - initial_time
    logger().info('Finished running simulation {} s'.format(elapsed_time))

    logger().info('Discarding initial transient equilibration...')
    # Discard initial transient equilibration
    [t0, g, Neff_max] = timeseries.detectEquilibration(u_n)
    x_n = x_n[t0:, :, :]
    u_n = u_n[t0:]

    # subsample to remove correlation
    logger().info('Subsample to remove correlation...')
    indices = timeseries.subsampleCorrelatedData(u_n, g=g)
    x_n = x_n[indices, :, :]
    u_n = u_n[indices]

    # store data
    if solvated:
        database['solvated']['context'] = context
        database['solvated']['x_n'] = x_n
        database['solvated']['u_n'] = u_n
    else:
        database['vacuum']['context'] = context
        database['vacuum']['x_n'] = x_n
        database['vacuum']['u_n'] = u_n

    # clean up
    del context, integrator

    logger().info("simulation %12.3f s | %5d samples discarded | %5d independent samples remain" % (elapsed_time, t0, len(indices)))

    return database


def _determine_fastest_platform():
    """
    Return the fastest available platform.
    Returns
    -------
    platform : simtk.openmm.Platform
    The fastest available platform.
    """
    platform_speeds = np.array([mm.Platform.getPlatform(i).getSpeed()
                                    for i in range(mm.Platform.getNumPlatforms())])
    fastest_platform_id = int(np.argmax(platform_speeds))
    platform = mm.Platform.getPlatform(fastest_platform_id)
    return platform


def zwanzig(database):

    temperature = 300*u.kelvin
    kT = kB * temperature
    beta = 1.0 / kT

    system = database['system']
    context = database['context']
    structure = database['structure']
    x_n =database['x_n']
    u_1 = database['u_n']
    n_samples = len(u_1)

    # copy torsions
    forces = {system.getForce(i).__class__.__name__: system.getForce(i)
                  for i in range(system.getNumForces())}
    torsion_force = forces['PeriodicTorsionForce']

    # create new force
    new_torsion_force = structure.omm_dihedral_force()
    # copy parameters
    for i in range(new_torsion_force.getNumTorsions()):
        torsion = new_torsion_force.getTorsionParameters(i)
        torsion_force.setTorsionParameters(i, *torsion)
    # update parameters in context
    torsion_force.updateParametersInContext(database['context'])
    #clean up
    del new_torsion_force

    # caclulate energy with new parameters for each sample
    u_2 = np.zeros(n_samples)
    for sample in range(n_samples):
        positions = u.Quantity(x_n[sample, :, :], u.nanometers)
        context.setPositions(positions)
        state = context.getState(getEnergy=True)
        u_2[sample] = beta * state.getPotentialEnergy()

    # calculate energy difference
    diff = u_1 - u_2
    expectation = np.average(np.exp(-diff))
    ddf_12 = kT * np.log(expectation)

    return ddf_12


def create_context(database, param, platform=None):
    #self.structure.load_parameters(param, copy_parameters=False)
    #self.system = self.structure.createSystem()
    platform = mm.Platform.getPlatformByName('Reference')
    time_step = 2.0*u.femtoseconds
    temperature = 300*u.kelvin
    friction_coef = 1.0/u.picoseconds
    integrator_solv = mm.LangevinIntegrator(temperature, friction_coef, time_step)
    integrator_vacuum = mm.LangevinIntegrator(temperature, friction_coef, time_step)


    solvated_system = database['solvated']['system']
    vacuum_system = database['vacuum']['system']
    if platform != None:
        database['solvated']['context'] = mm.Context(database['solvated']['system'], integrator_solv, platform)
        database['vacuum']['context'] = mm.Context(database['vacuum']['system'], integrator_vacuum, platform)


# def _opencl_device_support_precision(precision_model):
#     """
#     Check if this device supports the given precision model for OpenCL platform.
#     Some OpenCL devices do not support double precision. This offers a test
#     function.
#     Returns
#     -------
#     is_supported : bool
#     True if this device supports double precision for OpenCL, False
#     otherwise.
#     """
#     opencl_platform = mm.Platform.getPlatformByName('OpenCL')
#
#     # Platforms are singleton so we need to store
#     # the old precision model before modifying it
#     old_precision = opencl_platform.getPropertyDefaultValue('OpenCLPrecision')
#
#     # Test support by creating a toy context
#     YamlBuilder._set_gpu_precision(opencl_platform, precision_model)
#     system = mm.System()
#     system.addParticle(1.0 * u.amu)  # system needs at least 1 particle
#     integrator = mm.VerletIntegrator(1.0 * u.femtoseconds)
#     try:
#         context = mm.Context(system, integrator, opencl_platform)
#         is_supported = True
#     except Exception:
#         is_supported = False
#     else:
#         del context
#     del integrator
#
#     # Restore old precision
#     YamlBuilder._set_gpu_precision(opencl_platform, old_precision)
#
#     return is_supported




















