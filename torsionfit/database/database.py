# -*- coding: utf-8 -*-
"""
Created on Thu May  7 19:32:22 2015

"""

__author__ = 'Chaya D. Stern'

import numpy as np

import simtk.openmm as mm
from simtk.unit import Quantity, nanometers, kilojoules_per_mole, picoseconds

from  mdtraj import Trajectory
from torsionfit import parameters as par


class DataBase(Trajectory):
    """container object for molecular configurations and energies.

    Attributes
    ----------
    structure: ParmEd.Structure
    mm_energy: simtk.unit.Quantity((n_frames), unit=kilojoule/mole)
    positions:
    context:
    system:
    integrator:
    """

    def __init__(self, positions, topology, structure, time=None):
        """Create new TorsionScanSet object"""
        assert isinstance(topology, object)
        super(DataBase, self).__init__(positions, topology, time)
        self.structure = structure
        self.mm_energy = Quantity()
        self.positions = positions
        self.context = None
        self.system = None
        self.integrator = mm.VerletIntegrator(0.004*picoseconds)

        # Don't allow an empty TorsionScanSet to be created
        if self.n_frames == 0:
            msg = 'DataBase has no frames!\n'
            msg += '\n'
            msg += 'DataBase provided were:\n'
            msg += str(positions)
            raise Exception(msg)

    def create_context(self, param, platform=None):
        """

        Parameters
        ----------
        param :
        platform :
        """
        self.structure.load_parameters(param, copy_parameters=False)
        self.system = self.structure.createSystem()
        if platform != None:
            self.context = mm.Context(self.system, self.integrator, platform)
        else:
            self.context = mm.Context(self.system, self.integrator)

    def copy_torsions(self, param=None, platform=None):
        """

        Parameters
        ----------
        param :
        platform :
        """
        forces = {self.system.getForce(i).__class__.__name__: self.system.getForce(i)
                  for i in range(self.system.getNumForces())}
        torsion_force = forces['PeriodicTorsionForce']

        # create new force
        new_torsion_force = self.structure.omm_dihedral_force()

        # sanity check
        if torsion_force.getNumTorsions() != new_torsion_force.getNumTorsions():
            # create new context and new integrator. First delete old context and integrator
            del self.system
            del self.context
            del self.integrator
            self.integrator = mm.VerletIntegrator(0.004*u.picoseconds)
            self.create_context(param, platform)
            forces = {self.system.getForce(i).__class__.__name__: self.system.getForce(i)
                      for i in range(self.system.getNumForces())}
            torsion_force = forces['PeriodicTorsionForce']

        # copy parameters
        for i in range(new_torsion_force.getNumTorsions()):
            torsion = new_torsion_force.getTorsionParameters(i)
            torsion_force.setTorsionParameters(i, *torsion)
        # update parameters in context
        torsion_force.updateParametersInContext(self.context)

        # clean up
        del new_torsion_force

    def _string_summary_basic(self):
        """Basic summary of TorsionScanSet in string form."""
        energy_str = 'with MM Energy' if self._have_mm_energy else 'without MM Energy'
        value = "Database with %d frames, %d atoms, %d residues, %s" % (
                     self.n_frames, self.n_atoms, self.n_residues, energy_str)
        return value

    def energy(self, param, platform=None):
        """ Computes energy for a given structure with a given parameter set

        Parameters
        ----------
        offset :
        param: parmed.charmm.CharmmParameterSet
        platform: simtk.openmm.Platform to evaluate energy on (if None, will select automatically)
        """

        if self.n_frames == 0:
            raise Exception("self.n_frames = 0! There are no frames to compute energy for.")

        # Check if context exists.
        if not self.context:
            self.create_context(param, platform)
        else:
            # copy new torsion parameters
            self.copy_torsions(param, platform)

        # # Save initial mm energy
        # save = False
        # if not self._have_mm_energy:
        #     save = True
        # Compute potential energies for all snapshots.
        self.mm_energy = Quantity(value=np.zeros([self.n_frames], np.float64), unit=kilojoules_per_mole)
        for i in range(self.n_frames):
            self.context.setPositions(self.positions[i])
            state = self.context.getState(getEnergy=True)
            self.mm_energy[i] = state.getPotentialEnergy()

        # # Subtract off minimum of mm_energy and add offset
        # energy_unit = kilojoules_per_mole
        #
        # min_energy = self.mm_energy.min()
        # self.mm_energy -= min_energy
        # if save:
        #     self.initial_mm = deepcopy(self.mm_energy)
        # if offset:
        #     offset = Quantity(value=offset.value, unit=energy_unit)
        #     self.mm_energy += offset
        # self.delta_energy = (self.qm_energy - self.mm_energy)
        # # self.delta_energy = self.delta_energy - self.delta_energy.min()

    def mm_from_param_sample(self, param, db, start=0, end=-1, decouple_n=False, phase=False):

        """
        This function computes mm_energy for scan using sampled torsions
        Args:
            param: parmed.charmm.parameterset
            db: sqlit_plus database
            start: int, start of mcmc chain. Defualt 0
            end: int, end of mcmc chain. Default -1
            decouple_n: flag if multiplicities were sampled
            phase: flag if phases were sampled

        Returns:

        """
        N = len(db.sigma[start:end])
        mm_energy = np.zeros((N, self.n_frames))
        for i in range(N):
            par.param_from_db(param, db,  i, decouple_n, phase)
            self.compute_energy(param)
            mm_energy[i] = self.mm_energy._value

        return mm_energy

    @property
    def _have_mm_energy(self):
        return len(self.mm_energy) is not 0

    def __getitem__(self, key):
        "Get a slice of this trajectory"
        return self.slice(key)
