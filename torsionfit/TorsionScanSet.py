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
from cclib.parser import Gaussian, Psi
from cclib.parser.utils import convertor
from mdtraj import Trajectory
from simtk.unit import Quantity, nanometers, kilojoules_per_mole
from parmed.charmm import CharmmPsfFile, CharmmParameterSet
import torsionfit.utils as utils
from fnmatch import fnmatch
import os
import re


def to_optimize(param, stream, penalty=10):
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
    params_to_optimize = [k for k in param.dihedral_types.keys()
            if k not in keys and param.dihedral_types[k].penalty >= penalty]
    # compress type list
    written = set()
    for t in params_to_optimize:
            if t in written or tuple(reversed(t)) in written: continue
            written.add(t)
    return list(written)


def parse_psi4_log(logfiles, structure):
    """
    Parses output of psi4 torsion scan script
    :param logfiles: list of str
        logfiles of psi4 script
    :param structure: str
        Charmm psf file of structure
    :return:
    TorsionScanSet
    """
    topology = md.load_psf(structure)
    structure = CharmmPsfFile(structure)
    positions = np.ndarray((0, topology.n_atoms, 3))
    qm_energies = np.ndarray(0)
    torsions = np.ndarray((0, 4), dtype=int)
    directions = np.ndarray(0, dtype=int)
    angles = np.ndarray(0, dtype=float)

    if type(logfiles) != list:
        logfiles = [logfiles]

    for file in logfiles:
        qm = np.ndarray(0)
        fi = open(file, 'r')
        # check if log file is complete
        complete = False  # complete flag
        for line in fi:
            if line.startswith('Relative'):
                complete = True
        fi.seek(0)
        section = None
        torsion = np.ndarray((1, 4), dtype=int)
        angle = np.ndarray(0, dtype=float)
        for line in fi:
            # Flag if structure is optimized
            optimized = False
            if line.startswith('Optimizer'):
                optimized = True
                # Store Dihedral and position of optimized structures
                fi.next()
                l = filter(None, fi.next().strip().split(' '))
                dih = round(float(l[-2]))
                try:
                    t = l[-6:-2]
                    for i in range(len(t)):
                        torsion[0][i] = int(t[i]) - 1
                    torsions = np.append(torsions, torsion, axis=0)
                except ValueError:
                    pass
                angle = np.append(angle, dih)
                fi.next()
                pos = filter(None, re.split("[, \[\]]", fi.next().strip()))
                pos = [float(i) for i in pos]
                pos = np.asarray(pos).reshape((-1, 3))
                # convert angstroms to nanometers
                positions = np.append(positions, pos[np.newaxis]*0.1, axis=0)
            if not complete and optimized:
                # Find line that starts with energy
                for line in fi:
                    if line.startswith('Energy'):
                        energy = filter(None, line.strip().split(' '))[-1]
                        # Convert to KJ/mol
                        energy = float(energy)*2625.5
                        qm = np.append(qm, energy)
                        break
            if line.startswith('Relative'):
                section = 'Energy'
                fi.next()
                continue
            if section == 'Energy':
                line = filter(None, line.strip().split(' '))
                if line != []:
                    dih = round(float(line[0]))
                    if dih in angle:
                        # Only save energies of optimized structures
                        qm_energies = np.append(qm_energies, float(line[-1]))
        if qm.size is not 0:
            qm = qm - min(qm)
            qm_energies = np.append(qm_energies, qm)

        fi.close()
        angles = np.append(angles, angle, axis=0)
    return TorsionScanSet(positions, topology, structure, torsions, directions, angles, qm_energies)


def parse_psi4_out(oufiles_dir, structure):
    """
    Parse psi4 out files from distributed torsion scan (there are many output files, one for each structure)
    :param oufiles_dir: str
        path to directory where the psi4 output files are
    :param structure: str
        path to psf file of structure
    :return: TorsionScanSet

    """
    topology = md.load_psf(structure)
    structure = CharmmPsfFile(structure)
    positions = np.ndarray((0, topology.n_atoms, 3))
    qm_energies = np.ndarray(0)
    directions = np.ndarray(0, dtype=int)
    torsions = np.ndarray((0, 4), dtype=int)
    angles = np.ndarray(0, dtype=float)
    optimized = np.ndarray(0, dtype=bool)

    out_files = []
    dih_angle = []
    pattern = "*.out2"
    for path, subdir, files in os.walk(oufiles_dir):
        for name in files:
            if fnmatch(name, pattern):
                if name.startswith('timer'):
                    continue
                dih_angle.append(int(name.split('_')[-1].split('.')[0]))
                path = os.path.join(os.getcwd(), path, name)
                out_files.append(path)
    # Sort files in increasing angles order
    out_files = [out_file for (angle, out_file) in sorted(zip(dih_angle, out_files))]
    dih_angle.sort()
    if not out_files:
        raise Exception("There are no psi4 output files. Did you choose the right directory?")

    # Parse files
    for f in out_files:
        torsion = np.ndarray((1, 4), dtype=int)
        fi = open(f, 'r')
        for line in fi:
            if line.startswith('dih_string'):
                t = line.strip().split('"')[1].split(' ')[:4]
                for i in range(len(t)):
                    torsion[0][i] = int(t[i]) - 1
                torsions = np.append(torsions, torsion, axis=0)
        fi.close()
        optimizer = True
        log = Psi(f)
        data = log.parse()
        try:
            data.optdone
        except AttributeError:
            optimizer = False
            print("Warning: Optimizer failed for {}".format(f))
        optimized = np.append(optimized, optimizer)

        positions = np.append(positions, data.atomcoords[-1][np.newaxis]*0.1, axis=0)
        # Try MP2 energies. Otherwise take SCFenergies
        try:
            qm_energy = convertor(data.mpenergies[-1], "eV", "kJmol-1")
        except AttributeError:
            qm_energy = convertor(data.scfenergies[-1], "eV", "kJmol-1")
        qm_energies = np.append(qm_energies, qm_energy, axis=0)

    # Subtract lowest energy to find relative energies
    qm_energies = qm_energies - min(qm_energies)
    angles = np.asarray(dih_angle)
    return TorsionScanSet(positions, topology, structure, torsions, directions, angles, qm_energies, optimized)


def parse_gauss(logfiles, structure):
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

    for file in (logfiles):
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
    structure: ParmEd.Structure
    qm_energy: simtk.unit.Quantity((n_frames), unit=kilojoule/mole)
    mm_energy: simtk.unit.Quantity((n_frames), unit=kilojoule/mole)
    delta_energy: simtk.unit.Quantity((n_frames), unit=kilojoule/mole)
    torsion_index: {np.ndarray, shape(n_frames, 4)}
    step: {np.ndarray, shape(n_frame, 3)}
    direction: {np.ndarray, shape(n_frame)}. 0 = negative, 1 = positive
    """

    def __init__(self, positions, topology, structure, torsions, directions, steps, qm_energies, optimized=None):
        """Create new TorsionScanSet object"""
        assert isinstance(topology, object)
        super(TorsionScanSet, self).__init__(positions, topology)
        self.structure = structure
        self.qm_energy = Quantity(value=qm_energies, unit=kilojoules_per_mole)
        self.mm_energy = Quantity()
        self.initial_mm = Quantity()
        self.delta_energy = Quantity()
        self.torsion_index = torsions
        self.direction = directions
        self.steps = steps
        self.positions = positions
        self.context = None
        self.system = None
        self.integrator = mm.VerletIntegrator(0.004*u.picoseconds)
        self.energy = np.array
        self.optimized = optimized

        # Don't allow an empty TorsionScanSet to be created
        if self.n_frames == 0:
            msg = 'TorsionScanSet has no frames!\n'
            msg += '\n'
            msg += 'positions provided were:\n'
            msg += str(positions)
            raise Exception(msg)

    def create_context(self, param, platform=None):
        self.structure.load_parameters(param, copy_parameters=False)
        self.system = self.structure.createSystem()
        if platform != None:
            self.context = mm.Context(self.system, self.integrator, platform)
        else:
            self.context = mm.Context(self.system, self.integrator)

    def copy_torsions(self, param=None, platform=None):
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

    def to_dataframe(self, psi4=True):
        """ convert TorsionScanSet to pandas dataframe """

        data = []
        if psi4:
            for i in range(self.n_frames):
                if len(self.mm_energy) == self.n_frames and len(self.delta_energy) == self.n_frames:
                    try:
                        data.append((self.torsion_index[i], self.steps[i], self.qm_energy[i], self.mm_energy[i], self.delta_energy[i], self.optimized[i]))

                    except IndexError:
                        data.append((float('nan'), self.steps[i], self.qm_energy[i], self.mm_energy[i], self.delta_energy[i], float('nan')))

                else:
                    data.append((self.torsion_index[i], self.steps[i], self.qm_energy[i], float('nan'), float('nan'), self.optimized[i]))
            torsion_set = pd.DataFrame(data, columns=["Torsion", "Torsion angle", "QM energy KJ/mol", "MM energy KJ/mole",
                                                          "delat KJ/mol", "Optimized"])

        else:
            for i in range(self.n_frames):
                if len(self.mm_energy) == self.n_frames and len(self.delta_energy) == self.n_frames:
                    data.append((self.torsion_index[i], self.direction[i], self.steps[i], self.qm_energy[i],
                                 self.mm_energy[i], self.delta_energy[i]))
                else:
                    data.append((self.torsion_index[i], self.direction[i], self.steps[i], self.qm_energy[i],
                                 float('nan'), float('nan')))

            torsion_set = pd.DataFrame(data, columns=[ "torsion", "scan_direction", "step_point_total",
                                                       "QM_energy KJ/mol", "MM_energy KJ/mole", "delta KJ/mole"])

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

    def compute_energy(self, param, offset=None, platform=None,):
        """ Computes energy for a given structure with a given parameter set

        Parameters
        ----------
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

        # Save initial mm energy
        save = False
        if not self._have_mm_energy:
            save = True
        # Compute potential energies for all snapshots.
        self.mm_energy = Quantity(value=np.zeros([self.n_frames], np.float64), unit=kilojoules_per_mole)
        for i in range(self.n_frames):
            self.context.setPositions(self.positions[i])
            state = self.context.getState(getEnergy=True)
            self.mm_energy[i] = state.getPotentialEnergy()

        # Subtract off minimum of mm_energy and add offset
        energy_unit = kilojoules_per_mole

        min_energy = self.mm_energy.min()
        self.mm_energy -= min_energy
        if save:
            self.initial_mm = self.mm_energy
        if offset:
            offset = Quantity(value=offset.value, unit=energy_unit)
            self.mm_energy += offset
        self.delta_energy = (self.qm_energy - self.mm_energy)
        self.delta_energy = self.delta_energy - self.delta_energy.min()

        # Compute deviation between MM and QM energies with offset
        #self.delta_energy = mm_energy - self.qm_energy + Quantity(value=offset, unit=kilojoule_per_mole)

    def mm_from_param_sample(self, param, db, start=0, end=-1, decouple_n=False, phase=False):

        """
        This function computes mm_energy for scan using sampled torsions
        Args:
            param:
            db:
            start:
            end:
            decouple_n:
            phase:

        Returns:

        """
        N = len(db.sigma[start:end])
        mm_energy = np.zeros((N, self.n_frames))
        for i in range(N):
            utils.param_from_db(param, db,  i, decouple_n, phase)
            self.compute_energy(param)
            mm_energy[i] = self.mm_energy._value

        return mm_energy

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
