__author__ = 'Chaya D. Stern'

import pandas as pd
import numpy as np

from simtk.unit import Quantity, nanometers, kilojoules_per_mole

from cclib.parser import Gaussian, Psi
from cclib.parser.utils import convertor

import mdtraj as md
from parmed.charmm import CharmmPsfFile, CharmmParameterSet
from torsionfit.database import DataBase

from copy import deepcopy
from fnmatch import fnmatch
import os
import re
import warnings
import itertools


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
    return QMDataBase(positions, topology, structure, torsions, directions, angles, qm_energies)


def parse_psi4_out(oufiles_dir, structure, pattern="*.out"):
    """
    Parse psi4 out files from distributed torsion scan (there are many output files, one for each structure)
    :param oufiles_dir: str
        path to directory where the psi4 output files are
    :param structure: str
        path to psf file of structure
    :param pattern: str
        pattern for psi4 output file. Default is *.out
    :return: TorsionScanSet

    """
    topology = md.load_psf(structure)
    structure = CharmmPsfFile(structure)
    positions = np.ndarray((0, topology.n_atoms, 3))
    qm_energies = np.ndarray(0)
    torsions = np.ndarray((0, 4), dtype=int)
    angles = np.ndarray(0, dtype=float)
    optimized = np.ndarray(0, dtype=bool)

    out_files = {}
    for path, subdir, files in os.walk(oufiles_dir):
        for name in files:
            if fnmatch(name, pattern):
                if name.startswith('timer'):
                    continue
                name_split = name.split('_')
                try:
                    torsion_angle = (name_split[1] + '_' + name_split[2] + '_' + name_split[3] + '_' + name_split[4])
                except IndexError:
                    warnings.warn("Do you only have one torsion scan? The output files will be treated as one scan")
                    torsion_angle = 'only_one_scan'
                try:
                    out_files[torsion_angle]
                except KeyError:
                    out_files[torsion_angle] = []
                path = os.path.join(os.getcwd(), path, name)
                out_files[torsion_angle].append(path)
    # Sort files in increasing angles order for each torsion
    sorted_files = []
    dih_angles = []
    for tor in out_files:
        dih_angle = []
        for file in out_files[tor]:
            dih_angle.append(int(file.split('_')[-1].split('.')[0]))
        sorted_files.append([out_file for (angle, out_file) in sorted(zip(dih_angle, out_files[tor]))])
        dih_angle.sort()
        dih_angles.append(dih_angle)
    if not out_files:
        raise Exception("There are no psi4 output files. Did you choose the right directory?")

    # Parse files
    for f in itertools.chain.from_iterable(sorted_files):
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
            warnings.warn("Warning: Optimizer failed for {}".format(f))
        optimized = np.append(optimized, optimizer)

        positions = np.append(positions, data.atomcoords[-1][np.newaxis]*0.1, axis=0)
        # Try MP2 energies. Otherwise take SCFenergies
        try:
            qm_energy = convertor(data.mpenergies[-1], "eV", "kJmol-1")
        except AttributeError:
            try:
                qm_energy = convertor(np.array([data.scfenergies[-1]]), "eV", "kJmol-1")
            except AttributeError:
                warnings.warn("Warning: Check if the file terminated before completing SCF")
                qm_energy = np.array([np.nan])
        qm_energies = np.append(qm_energies, qm_energy, axis=0)

    # Subtract lowest energy to find relative energies
    qm_energies = qm_energies - min(qm_energies)
    angles = np.asarray(list(itertools.chain.from_iterable(dih_angles)))
    return QMDataBase(positions=positions, topology=topology, structure=structure, torsions=torsions, angles=angles,
                      qm_energies=qm_energies, optimized=optimized)


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
    return QMDataBase(positions=positions, topology=topology, structure=structure, torsions=torsions, steps=steps,
                      qm_energies=qm_energies, directions=directions)


class QMDataBase(DataBase):
    """container object for torsion scan

    A TorsionScanSet should be constructed by loading Gaussian 09 torsion scan log files or a psi4 output file from disk
    with an mdtraj.Topology object

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

    def __init__(self, positions, topology, structure, torsions, qm_energies, angles=None, steps=None, directions=None,
                 optimized=None, time=None):
        """Create new TorsionScanSet object"""
        assert isinstance(topology, object)
        super(QMDataBase, self).__init__(positions, topology, structure, time)
        self.qm_energy = Quantity(value=qm_energies, unit=kilojoules_per_mole)
        self.initial_mm = Quantity()
        self.delta_energy = Quantity()
        self.torsion_index = torsions
        self.direction = directions
        self.steps = steps
        self.angles = angles
        self.optimized = optimized
        self.phis = {}

    def compute_energy(self, param, offset=None, platform=None):
        """ Computes energy for a given structure with a given parameter set

        Parameters
        ----------
        param: parmed.charmm.CharmmParameterSet
        platform: simtk.openmm.Platform to evaluate energy on (if None, will select automatically)
        """

        # Save initial mm energy
        save = False
        if not self._have_mm_energy:
            save = True

        # calculate energy
        super(QMDataBase, self).compute_energy(param, platform)

        # Subtract off minimum of mm_energy and add offset
        energy_unit = kilojoules_per_mole

        min_energy = self.mm_energy.min()
        self.mm_energy -= min_energy
        if save:
            self.initial_mm = deepcopy(self.mm_energy)
        if offset:
            offset = Quantity(value=offset.value, unit=energy_unit)
            self.mm_energy += offset
        self.delta_energy = (self.qm_energy - self.mm_energy)
        # self.delta_energy = self.delta_energy - self.delta_energy.min()

    def to_dataframe(self, psi4=True):

        """ convert TorsionScanSet to pandas dataframe

        Parameters
        ----------
        psi4 : bool
            Flag if QM log file is from psi4. Default True.
        """

        if len(self.mm_energy) == self.n_frames and len(self.delta_energy) == self.n_frames:
            mm_energy = self.mm_energy
            delta_energy = self.delta_energy
        else:
            mm_energy = [float('nan') for _ in range(self.n_frames)]
            delta_energy = [float('nan') for _ in range(self.n_frames)]
        if psi4:
            data = [(self.torsion_index[i], self.angles[i], self.qm_energy[i], mm_energy[i], delta_energy[i],
                    self.optimized[i]) for i in range(self.n_frames)]

            columns = ['Torsion', 'Torsion angle', 'QM energy (KJ/mol)', 'MM energy (KJ/mol)', 'Delta energy (KJ/mol)',
                       'Optimized']
        else:
            data = [(self.torsion_index[i], self.direction[i], self.steps[i], self.qm_energy[i], mm_energy[i],
                     delta_energy[i]) for i in range(self.n_frames)]
            columns = ['Torsion', 'direction','steps', 'QM energy (KJ/mol)', 'MM energy (KJ/mol)',
                       'Delta energy (KJ/mol)']

        torsion_set = pd.DataFrame(data, columns=columns)
        return torsion_set

    def extract_geom_opt(self):
        """
        Extracts optimized geometry for Gaussian torsion scan.

        Returns
        -------
        New QMDataBase

        """
        key = []
        for i, step in enumerate(self.steps):
            try:
                if step[1] != self.steps[i+1][1]:
                    key.append(i)
            except IndexError:
                key.append(i)
        new_torsionScanSet = self.slice(key)
        return new_torsionScanSet

    def remove_nonoptimized(self):
        """
        Remove configurations where optimizer failed
        Returns: copy of scan set with only optimized structures

        Returns
        -------
        new QMDataBase

        """
        key = []
        for i, optimized in enumerate(self.optimized):
            if optimized:
                key.append(i)
        new_torsionscanset = self.slice(key)
        return new_torsionscanset

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
        if self.direction is not None:
            direction = self.direction[key]
        if self.optimized is not None:
            optimized = self.optimized[key]
        if self.steps is not None:
            steps = self.steps[key]
        if self.angles is not None:
            angles = self.angles[key]
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
            qm_energy = qm_energy.copy()
            if self.direction is not None:
                direction = direction.copy()
            else:
                direction = self.direction
            if self.optimized is not None:
                    optimized = optimized.copy()
            else:
                optimized = self.optimized
            if self.steps is not None:
                steps = steps.copy()
            else:
                steps = self.steps
            if self.angles is not None:
                angles = angles.copy()
            else:
                angles = self.angles
            if self.unitcell_angles is not None:
                unitcell_angles = unitcell_angles.copy()
            if self.unitcell_lengths is not None:
                unitcell_lengths = unitcell_lengths.copy()

        newtraj = self.__class__(
            positions=xyz, topology=topology, structure=structure, torsions=torsions, directions=direction, steps=steps,
            qm_energies=qm_energy, optimized=optimized, angles=angles, time=time)

        if self._rmsd_traces is not None:
            newtraj._rmsd_traces = np.array(self._rmsd_traces[key],
                                            ndmin=1, copy=True)
        return newtraj

    # def combine(self, qmdabase, copy=True):
    #     """
    #     Add more QM configurations to database
    #
    #     Parameters
    #     ----------
    #     directory : path to directory of QM scan
    #
    #     Returns
    #     -------
    #     Updates database in place
    #
    #     """
    #
    #     xyz = self.xyz
    #     time = self.time
    #     torsions = self.torsion_index
    #     if self.direction is not None:
    #         direction = self.direction
    #     if self.optimized is not None:
    #         optimized = self.optimized
    #     if self.steps is not None:
    #         steps = self.steps
    #     if self.angles is not None:
    #         angles = self.angles
    #     qm_energy = self.qm_energy
    #     unitcell_lengths, unitcell_angles = None, None
    #     if self.unitcell_angles is not None:
    #         unitcell_angles = self.unitcell_angles
    #     if self.unitcell_lengths is not None:
    #         unitcell_lengths = self.unitcell_lengths
    #
    #     if copy:
    #         xyz = xyz.copy()
    #         time = time.copy()
    #         topology = deepcopy(self._topology)
    #         structure = deepcopy(self.structure)
    #         torsions = torsions.copy()
    #         qm_energy = qm_energy.copy()
    #         if self.direction is not None:
    #             direction = direction.copy()
    #         else:
    #             direction = self.direction
    #         if self.optimized is not None:
    #                 optimized = optimized.copy()
    #         else:
    #             optimized = self.optimized
    #         if self.steps is not None:
    #             steps = steps.copy()
    #         else:
    #             steps = self.steps
    #         if self.angles is not None:
    #             angles = angles.copy()
    #         else:
    #             angles = self.angles
    #         if self.unitcell_angles is not None:
    #             unitcell_angles = unitcell_angles.copy()
    #         if self.unitcell_lengths is not None:
    #             unitcell_lengths = unitcell_lengths.copy()
    #
    #     newtraj = self.__class__(
    #         positions=xyz, topology=topology, structure=structure, torsions=torsions, directions=direction, steps=steps,
    #         qm_energies=qm_energy, optimized=optimized, angles=angles, time=time)
    #
    #     if self._rmsd_traces is not None:
    #         newtraj._rmsd_traces = np.array(self._rmsd_traces[key],
    #                                         ndmin=1, copy=True)
    #     return newtraj

    def build_phis(self, to_optimize=None):
        """
        This function builds a dictionary of phis for specified dihedrals in the molecules for all frames in the qm db.

        Parameters
        ----------
        to_optimize : list of dihedral types to calculate phis for
            Default is None. When None, it will calculate phis for all dihedral types in molecule

        """

        type_list = to_optimize
        if type_list is None:
            type_list = []
            for torsion_type in self.structure.dihedrals:
                t = (torsion_type.atom1.type, torsion_type.atom2.type, torsion_type.atom3.type, torsion_type.atom4.type)
                type_list.append(t)

        type_frequency = {}
        for t in type_list:
            # If t is not a palindrome, reverse it.
            if t[0] >= t[-1]:
                t = tuple(reversed(t))
            try:
                type_frequency[t] += 1
            except KeyError:
                type_frequency[t] = 1

        # sanity check
        if len(self.structure.dihedrals) != sum(type_frequency.values()):
            warnings.warn("type frequency values don't sum up to number of dihedral")

        self.phis = {t_type: [[] for i in range(self.n_frames)] for t_type in type_frequency}
        for i in range(self.n_frames):
            for dihedral in self.structure.dihedrals:
                atom = dihedral.atom1
                bond_atom = dihedral.atom2
                angle_atom = dihedral.atom3
                torsion_atom = dihedral.atom4

                torsion_type = (atom.type, bond_atom.type, angle_atom.type, torsion_atom.type)
                try:
                    self._append_phi(i, torsion_type, atom, bond_atom, angle_atom, torsion_atom)
                except KeyError:
                    warnings.warn("torsion {} is not in list of phis to precalculate but is in the structure. "
                                  "Are you sure you did not want to fit it?".format(torsion_type))
                # try:
                #     self.phis[torsion_type][i].append(self._cartesian_to_phi(atom, bond_atom, angle_atom, torsion_atom,
                #                                                              i))
                # except KeyError:
                #     self.phis[tuple(reversed(torsion_type))][i].append(self._cartesian_to_phi(atom, bond_atom,
                #                                                        angle_atom, torsion_atom, i))

        # Convert to np.array
        for t in self.phis:
            self.phis[t] = np.array(self.phis[t])

    def _append_phi(self, i, torsion_type, atom, bond_atom, angle_atom, torsion_atom):
        """
        Helper function to try to append a calculated phi angle to existing list for a torsion type

        Parameters
        ----------
        i : int
            frame for which phi is being calculated.
        torsion_type : tuple of strings
        atom : parmed atom type
            first atom in dihedral
        bond_atom : parmed atomtype
            second atom in dihedral
        angle_atom : parmed atomtype
            third atom in dihedral
        torsion_atom : parmed atomtype
            fourth atom in dihedral

        """
        try:
            self.phis[torsion_type][i].append(self._cartesian_to_phi(atom, bond_atom, angle_atom, torsion_atom,
                                                                             i))
        except KeyError:
            self.phis[tuple(reversed(torsion_type))][i].append(self._cartesian_to_phi(atom, bond_atom,
                                                                       angle_atom, torsion_atom, i))

    def _cartesian_to_phi(self, atom, bond_atom, angle_atom, torsion_atom, i):
        """
        measures torsion angle for a specific torsion

        Parameters
        ----------
        atom : parmed atom
        bond_atom : parmed atom
        angle_atom : parmed atom
        torsion_atom : parmed atom
        i : int
            index for configuration (frame)

        Returns
        -------
        phi: float;
            torsion angle in radians

        """

        atom1_coords = self.positions[i][atom.idx]
        bond_coords = self.positions[i][bond_atom.idx]
        angle_coords = self.positions[i][angle_atom.idx]
        torsion_coords = self.positions[i][torsion_atom.idx]

        a = atom1_coords - bond_coords
        b = angle_coords - bond_coords
        #3-4 bond
        c = angle_coords - torsion_coords
        a_u = a / np.linalg.norm(a)
        b_u = b / np.linalg.norm(b)
        c_u = c / np.linalg.norm(c)

        plane1 = np.cross(a_u, b_u)
        plane2 = np.cross(b_u, c_u)

        cos_phi = np.dot(plane1, plane2) / (np.linalg.norm(plane1)*np.linalg.norm(plane2))
        cos_phi = np.dot(plane1, plane2) / (np.linalg.norm(plane1)*np.linalg.norm(plane2))
        if cos_phi < -1.0:
            cos_phi = -1.0
        elif cos_phi > 1.0:
            cos_phi = 1.0

        phi = np.arccos(cos_phi)

        if np.dot(a, plane2) <= 0:
            phi = -phi

        return phi


    #def Fourier_series(self, K, n, phase):