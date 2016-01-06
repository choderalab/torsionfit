# -*- coding: utf-8 -*-
"""
Created on Thu May  7 19:32:22 2015

@author: sternc1
"""
import pandas as pd
import numpy as np
import simtk.openmm as mm
import simtk.unit as u
import mdtraj as md
from copy import copy, deepcopy
import re
from cclib.parser import Gaussian
from cclib.parser.utils import convertor 
from mdtraj import Trajectory
from simtk.unit import Quantity, nanometers, kilojoules_per_mole
from parmed.charmm import CharmmPsfFile


def to_optimize(param, stream, penalty = 10):
    """ returns a list of dihedrals to optimize and updates CharmmParameterSet
    with stream files

    Parameters
    ----------
    param : CharmmParameterSet
    stream: list of stream files
    penalty: int for CGenFF penalty cutoff (Default = 10)

    Returns list of tuples containing dihedrals to optimize

    """
    if type(stream) != list:
        stream = [stream]
    keys = [i for i in param.dihedral_types.keys()]
    for j in stream:
        param.read_stream_file(j)
    return [k for k in param.dihedral_types.keys()
            if k not in keys and param.dihedral_types[k].penalty >= penalty]


def read_scan_logfile(logfiles, structure):
    """ parses Guassian09 torsion-scan log file

    parameters
    ----------
    logfiles: str of list of str
                Name of Guassian 09 torsion scan log file
    structure: charmm psf file

    returns
    -------
    TorsionScanSet
    """
    topology = md.load_psf(structure)
    structure = CharmmPsfFile(structure)
    positions = np.ndarray((0, topology.n_atoms, 3))
    qm_energies = np.ndarray(0)
    torsions = np.ndarray((0, 4), dtype=int)
    directions = np.ndarray(0, dtype=int)
    steps = np.ndarray((0, 3), dtype=int)

    if type(logfiles) != list:
        logfiles = [logfiles]

    for file in sorted(logfiles):
        #print("loading %s" % file)
        direction = np.ndarray(1)
        torsion = np.ndarray((1, 4), dtype=int)
        step = np.ndarray((0, 3), dtype=int)
        index = (2, 12, -1)
        # f = file.split('/')[-1].split('.')
        # if f[2] == 'pos':
        #     direction[0] = 1
        # else:
        #     direction[0] = 0


        log = Gaussian(file)
        data = log.parse()
        # convert angstroms to nanometers
        positions = np.append(positions, data.atomcoords*0.1, axis=0)
        # Only add qm energies for structures that converged (because cclib throws out those coords but not other info)
        qm_energies = np.append(qm_energies, (convertor(data.scfenergies[:len(data.atomcoords)], "eV", "kJmol-1") -
                                              min(convertor(data.scfenergies[:len(data.atomcoords)], "eV", "kJmol-1"))), axis=0)

        fi = open(file, 'r')
        for line in fi:
            if re.search('   Scan   ', line):
                t = line.split()[2].split(',')
                t[0] = t[0][-1]
                t[-1] = t[-1][0]
                for i in range(len(t)):
                    torsion[0][i] = (int(t[i]) - 1)
            if re.search('^ D ', line):
                d = line.split()[-1]
                if d[0] == '-':
                    direction[0] = 0
                elif d[0] == '1':
                    direction[0] = 1
            if re.search('Step', line):
                try:
                    point = np.array(([int(line.rsplit()[j]) for j in index]))
                    point = point[np.newaxis,:]
                    step = np.append(step, point, axis=0)
                except:
                    pass

        fi.close()
        # only add scan points from converged structures
        steps = np.append(steps, step[:len(data.atomcoords)], axis=0)
        for i in range(len(data.atomcoords)):
            torsions = np.append(torsions, torsion, axis=0)
            directions = np.append(directions, direction, axis=0)
        del log
        del data
    return TorsionScanSet(positions, topology, structure, torsions, directions, steps, qm_energies)


class TorsionScanSet(Trajectory):
    """container object for torsion scan

    A TorsionScanSet should be constructed by loading Gaussian 09 torsion scan log files from disk
    with an mdtraj.Topology object

    Examples
    --------
    >>> torsion_set = read_scan_logfile('../examples/data/pyrrole/torsion-scan/PRL.scan2.neg.log', '../examples/data/pyrrole/pyrrol.psf')
    >>> print(torsion_set)
    <torsions.TorsionScanSet with 40 frames, 22 atoms, 1 residues, without MM Energy>


    Attributes
    ----------
    structure: chemistry.Structure
    qm_energy: simtk.unit.Quantity((n_frames), unit=kilojoule/mole)
    mm_energy: simtk.unit.Quantity((n_frames), unit=kilojoule/mole)
    delta_energy: simtk.unit.Quantity((n_frames), unit=kilojoule/mole)
    torsion_index: {np.ndarray, shape(n_frames, 4)}
    step: {np.ndarray, shape(n_frame, 3)}
    direction: {np.ndarray, shape(n_frame)}. 0 = negative, 1 = positive
    """

    def __init__(self, positions, topology, structure, torsions, directions, steps, qm_energies):
        """Create new TorsionScanSet object"""
        assert isinstance(topology, object)
        super(TorsionScanSet, self).__init__(positions, topology)
        self.structure = structure
        self.qm_energy = Quantity(value=qm_energies, unit=kilojoules_per_mole)
        self.mm_energy = Quantity()
        self.delta_energy = Quantity()
        self.torsion_index = torsions
        self.direction = directions
        self.steps = steps


    # def create_omm_system(self, param):
    #     """ Creates an OpenMM system for a given param set
    #     :param param: parmed.charmm.CharmmParameterSet
    #     :return: TorsionScanSet with omm system for each configuration
    #     """
    #
    #     # Create system
    #     system = self.structure.createSystem(param)
    #     self.system = system


    def to_dataframe(self):
        """ convert TorsionScanSet to pandas dataframe """

        data = []
        for i in range(self.n_frames):
            if len(self.mm_energy) == self.n_frames and len(self.delta_energy) == self.n_frames:
                data.append((self.torsion_index[i], self.direction[i], self.steps[i], self.qm_energy[i], self.mm_energy[i],
                             self.delta_energy[i]))
            else:
                data.append((self.torsion_index[i], self.direction[i], self.steps[i], self.qm_energy[i], float('nan'), float('nan')))

        torsion_set = pd.DataFrame(data, columns=[ "torsion", "scan_direction", "step_point_total", "QM_energy KJ/mol",
                                                   "MM_energy KJ/mole", "delta KJ/mole"])

        return torsion_set

    def _string_summary_basic(self):
        """Basic summary of TorsionScanSet in string form."""
        energy_str = 'with MM Energy' if self._have_mm_energy else 'without MM Energy'
        value = "torsions.TorsionScanSet with %d frames, %d atoms, %d residues, %s" % (
                     self.n_frames, self.n_atoms, self.n_residues, energy_str)
        return value

    def extract_geom_opt(self):
        key = []
        for i, step in enumerate(self.steps):
            try:
                if step[1] != self.steps[i+1][1]:
                    key.append(i)
            except IndexError:
                key.append(i)
        new_torsionScanSet = self.slice(key)
        return new_torsionScanSet

    def compute_energy(self,param, offset, platform=None,):
        """ Computes energy for a given structure with a given parameter set

        Parameters
        ----------
        param: chemistry.charmm.CharmmParameterSet
        platform: simtk.openmm.Platform to evaluate energy on (if None, will select automatically)
        """
        # Create Context.
        integrator = mm.VerletIntegrator(0.004*u.picoseconds)
        system = self.structure.createSystem(param)
        if platform != None:
            context = mm.Context(system, integrator, platform)
        else:
            context = mm.Context(system, integrator)

        # Compute potential energies for all snapshots.
        self.mm_energy = Quantity(value=np.zeros([self.n_frames], np.float64), unit=kilojoules_per_mole)
        for i in range(self.n_frames):
            context.setPositions(self.openmm_positions(i))
            state = context.getState(getEnergy=True)
            self.mm_energy[i] = state.getPotentialEnergy()

        # Subtract off minimum of mm_energy
        self.mm_energy -= self.mm_energy.min() + Quantity(value=float(offset.value), unit=kilojoules_per_mole)
        self.delta_energy = (self.qm_energy - self.mm_energy)

        # Compute deviation between MM and QM energies with offset
        #self.delta_energy = mm_energy - self.qm_energy + Quantity(value=offset, unit=kilojoule_per_mole)

        # Clean up.
        del context
        del system
        del integrator
        # print('Heap at end of compute_energy'), hp.heeap()

    @property
    def _have_mm_energy(self):
        return len(self.mm_energy) is not 0

    # @property
    # def _unique_torsions(self):
    # Not returning the right amount. debug
    #     torsions = []
    #     for i in range(len(self.torsion_index)):
    #         try:
    #             if (self.torsion_index[i] != self.torsion_index[i+1]).all():
    #                 torsions.append(self.torsion_index[i]), torsions.append(self.torsion_index[i+1])
    #         except:
    #             pass
    #     return len(torsions), torsions


    def __getitem__(self, key):
        "Get a slice of this trajectory"
        return self.slice(key)

    def slice(self, key, copy=True):
        """Slice trajectory, by extracting one or more frames into a separate object

        This method can also be called using index bracket notation, i.e
        `traj[1] == traj.slice(1)`

        Parameters
        ----------
        key : {int, np.ndarray, slice}
            The slice to take. Can be either an int, a list of ints, or a slice
            object.
        copy : bool, default=True
            Copy the arrays after slicing. If you set this to false, then if
            you modify a slice, you'll modify the original array since they
            point to the same data.
        """
        xyz = self.xyz[key]
        time = self.time[key]
        torsions = self.torsion_index[key]
        direction = self.direction[key]
        steps = self.steps[key]
        qm_energy = self.qm_energy[key]
        unitcell_lengths, unitcell_angles = None, None
        if self.unitcell_angles is not None:
            unitcell_angles = self.unitcell_angles[key]
        if self.unitcell_lengths is not None:
            unitcell_lengths = self.unitcell_lengths[key]

        if copy:
            xyz = xyz.copy()
            time = time.copy()
            topology = deepcopy(self._topology)
            structure = deepcopy(self.structure)
            torsions = torsions.copy()
            direction = direction.copy()
            steps = steps.copy()
            qm_energy = qm_energy.copy()

            if self.unitcell_angles is not None:
                unitcell_angles = unitcell_angles.copy()
            if self.unitcell_lengths is not None:
                unitcell_lengths = unitcell_lengths.copy()

        newtraj = self.__class__(
            xyz, topology, structure, torsions, direction, steps, qm_energy)

        if self._rmsd_traces is not None:
            newtraj._rmsd_traces = np.array(self._rmsd_traces[key],
                                            ndmin=1, copy=True)
        return newtraj




