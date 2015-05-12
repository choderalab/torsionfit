# -*- coding: utf-8 -*-
"""
Created on Thu May  7 19:32:22 2015

@author: sternc1
"""
from collections import OrderedDict
import pandas as pd
from cclib.parser import Gaussian
from cclib.parser.utils import convertor 
import re
from mdtraj import Trajectory
import numpy as np
import simtk.openmm as mm
from simtk.unit import Quantity, kilojoules_per_mole, nanometer


def read_scan_logfile(logfiles, topology):
    """
    parameters
    ----------
    logfiles: str of list of str
                Name of Guassian 09 torsion scan log file
    topology: mdtraj.Topology

    returns
    -------
    TorsionScanSet
    """

    positions = np.ndarray((0, topology.n_atoms, 3))
    energies = np.ndarray(0)
    torsions = np.ndarray((0, 4), dtype=int)
    directions = np.ndarray((0, 1), dtype=int)
    steps = np.ndarray((0, 3), dtype=int)

    if type(logfiles) != list:
        logfiles = [logfiles]

    for file in logfiles:
        print file
        direction = np.ndarray((1,1))
        torsion = np.ndarray((1,4), dtype=int)
        step = []
        index = (2, 12, -1)
        f = file.split('/')[-1].split('.')
        print f[2]
        if f[2] == 'pos':
            direction[0][0] = 1
        else:
            direction[0][0] = 0

        fi = open(file, 'r')
        for line in fi:

            if re.search('   Scan   ', line):
                t = line.split()[2].split(',')
                t[0] = t[0][-1]
                t[-1] = t[-1][0]
                for i in range(len(t)):
                    torsion[0][i] = (int(t[i]) - 1)
            if re.search('Step', line):
                try:
                    step = np.array(([int(line.rsplit()[j]) for j in index]))
                    step = step[np.newaxis,:]
                    steps = np.append(steps, step, axis=0)
                except:
                    pass
        fi.close()

        log = Gaussian(file)
        data = log.parse()
        # convert angstroms to nanometers
        positions = np.append(positions, data.atomcoords*0.1, axis=0)
        energies = np.append(energies, convertor(data.scfenergies, "eV", "kJmol-1"), axis=0)
        for i in range(len(data.scfenergies)):
            torsions = np.append(torsions, torsion, axis=0)
            directions = np.append(directions, direction, axis=0)

    qm_energy = Quantity(value=energies, unit=kilojoules_per_mole)

    return TorsionScanSet(positions, topology, torsions, directions, steps, qm_energy)


class TorsionScanSet(Trajectory):
    """container object for torsion scan
    
    A TorsionScanSet should be constructed by loading Gaussian 09 torsion scan log files from disk
    with an mdtraj.Topology object

    Examples
    --------
    >>> torsion_set = read_scan_logfile('FRG.scanN.dir.log')
    >>> print torsion_set
    <torsions.TorsionScanSet with 346 frames, 22 atoms, 1 residues, 4 unique torsions without MM Energy at 0x10b099b10>


    Attributes
    ----------
    qm_energy: simtk.unit.Quantity((n_frames, 1), unit=kilojoule/mole)
    mm_energy: simtk.unit.Quantity((n_frames, 1), unit=kilojoule/mole)
    delta_energy: simtk.unit.Quantity((n_frames, 1), unit=kilojoule/mole)
    torsion_index: {np.ndarray, shape(n_frames, 4)}
    step: {np.ndarray, shape(n_frame, 3)}
    direction: {np.ndarray, shape(n_frame, 1)}
    """

    def __init__(self, positions, topology, torsions, directions, steps, qm_energy):
        """Create new TorsionScanSet object"""
        Trajectory.__init__(self, positions, topology)
        self.qm_energy = qm_energy
        self.mm_energy = Quantity()
        self.delta_energy = Quantity()
        self.torsion_index = torsions
        self.direction = directions
        self.steps = steps

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
        """Basic summary of traj in string form."""
        energy_str = 'and MM Energy' if self._have_mm_energy else 'without MM Energy'
        value = "torsions.TorsionScanSet with %d frames, %d atoms, %d residues, %d unique torsions %s" % (
                     self.n_frames, self.n_atoms, self.n_residues, self._unique_torsions[0], energy_str)
        return value

    @property
    def _have_mm_energy(self):
        return len(self.mm_energy) is not 0

    @property
    def _unique_torsions(self):
        torsions = []
        for i in range(len(self.torsion_index)):
            try:
                if (self.torsion_index[i] != self.torsion_index[i+1]).all():
                    torsions.append(self.torsion_index[i]), torsions.append(self.torsion_index[i+1])
            except:
                pass
        return len(torsions), torsions





