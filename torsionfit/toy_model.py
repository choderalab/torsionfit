"""
Toy model of 4 carbons to test torsionfit

"""

from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed.topologyobjects import DihedralType, DihedralTypeList, Dihedral
import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as units
import mdtraj as md
import numpy as np
import torsionfit.TorsionScanSet as ScanSet
import torsionfit.TorsionFitModel as Model
from torsionfit.tests.utils import get_fun


class ToyModel(object):
    """
    Toy model of 4 carbons to test torsionfit

    Attributes
    ----------
    synthetic_energy: simtk.units.Quantity((n_increments), unit=kilojoule/mole)
    true_value: parmed.topologyobjects.DihedralTypeList (the dihedral parameters used to calculate synthetic energy)
    initital_value: parmed.topologyobjects.DihedralTypeList (dihedral parameters used to initialize pymc model)
    torsion_scan: torsionfit.TorsionScanSet
    model: torsionfit.TorsionFitModel

    The default toy model has randomized true and initial dihedral parameters. This can be changed by passing a
    parmed.topologyobjects.DihedralTypeList to true_value and/or initial_value when initializing a model instance.

    """

    def __init__(self, true_value='random', initial_value='random', phase='symmetric', multiplicity='on',
                 n_increments=13):
        self._param = CharmmParameterSet(get_fun('toy.str'))
        self._struct = CharmmPsfFile(get_fun('toy.psf'))
        self._pdb = app.PDBFile(get_fun('toy.pdb'))
        self._topology = md.load_psf(get_fun('toy.psf'))
        self.synthetic_energy = units.Quantity()
        self._positions = units.Quantity()
        self._platform = mm.Platform.getPlatformByName('Reference')

        # Replace ('CG331', 'CG321', 'CG321', 'CG331') torsion with true_value
        self._dih_type = ('CG331', 'CG321', 'CG321', 'CG331')
        original_torsion = self._param.dihedral_types[self._dih_type]
        if true_value == 'random':
            self._randomize_dih_param()
            self.true_value = self._param.dihedral_types[self._dih_type]
        else:
            dih_tlist = DihedralTypeList()
            for dih in true_value:
                dih_tlist.append(DihedralType(dih))
            self._param.dihedral_types[self._dih_type] = dih_tlist
            self.true_value = self._param.dihedral_types[self._dih_type]

        # parametrize toy
        self._struct.load_parameters(self._param, copy_parameters=False)
        self._struct.positions = self._pdb.positions

        # generate synthetic torsion scan
        self._torsion_scan(n_increments=n_increments)

        # initialize parameter
        if initial_value == 'random':
            self._randomize_dih_param()
            self.initial_value = self._param.dihedral_types[self._dih_type]
        else:
            # return original dihedral to param
            self._param.dihedral_types[self._dih_type] = original_torsion
            self.initial_value = self._param.dihedral_types[self._dih_type]

        # create torsionfit.TorsionScanSet
        torsions = np.zeros((len(self._positions), 4))
        torsions[:] = [1,2,3,4]
        direction = None
        steps = None
        self.scan_set = ScanSet.TorsionScanSet(self._positions.value_in_unit(units.nanometers), self._topology, self._struct,
                                  torsions, direction, steps, self.synthetic_energy.value_in_unit(units.kilojoules_per_mole))

        # create torsionfit.TorsionFitModel
        self.model = Model.TorsionFitModel(self._param, self._struct, self.scan_set, platform=self._platform,
                                           param_to_opt=[self._dih_type])

    def _spher2cart(self, r, theta, phi):
        """convert spherical to cartesian coordinates

        Paramters:
        r: bond length
        theta: bond angle
        phi: dihedral angle

        returns:
        cartesian coordinates
        """
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)
        return [x, y, z]

    def _torsion_scan(self, n_increments):
        """
        Generate positions and energy for torsion scan
        Parameters:
            n_increments: int. how many points to scan
        """
        n_increments = n_increments
        n_atoms = 4
        phis = np.arange(-np.pi/2, +np.pi/2, (np.pi)/n_increments)
        positions = np.zeros((len(phis), n_atoms, 3))

        # Get bond length in nm
        for bond in self._struct.bonds:
            if bond.atom1.type == bond.atom2.type:
                length = units.Quantity(value=bond.type.req * 0.1, unit=units.nanometers) # convert to nm
        # Get angle in radians
        for angle in self._struct.angles:
            if angle.atom1.type == angle.atom2.type:
                theta = units.Quantity(value=angle.type.theteq * (np.pi/180.0), unit=units.radians)

        atom1_coords = self._spher2cart(length.value_in_unit(units.nanometer), theta.value_in_unit(units.radian), phis[0])
        for i, phi in enumerate(phis):
            atom3_coords = self._spher2cart(length.value_in_unit(units.nanometer), theta.value_in_unit(units.radian), phi)
            atom3_coords[-1] = abs(atom3_coords[-1]) + length._value
            positions[i] = [atom1_coords,
                [0.000, 0.000, 0.000],
                [0.000, 0.000, length._value],
                atom3_coords]
        self._positions = units.Quantity(value=positions, unit=units.nanometer) # put units back in

        # calculate energy
        self.synthetic_energy = units.Quantity(value=np.zeros((len(positions))), unit=units.kilojoules_per_mole)
        platform = mm.Platform.getPlatformByName('Reference')
        integrator = mm.VerletIntegrator(0.004*units.picoseconds)
        system = self._struct.createSystem()
        context = mm.Context(system, integrator, platform)
        for i, conf in enumerate(self._positions):
            context.setPositions(conf)
            state = context.getState(getEnergy=True)
            self.synthetic_energy[i] = state.getPotentialEnergy()

    def _randomize_dih_param(self):
        """
        generates random dihedral parameters

        """
        dih_tlist = DihedralTypeList()
        multiplicities = [1, 2, 3, 4, 6]
        terms = np.random.randint(1, 5+1)
        np.random.shuffle(multiplicities)
        for i in range(terms):
            k = np.random.uniform(0.0, 20.0)
            n = multiplicities[i]
            phase = np.random.randint(0, 1+1)
            if phase == 1:
                phase = 180
            dih_tlist.append(DihedralType(k, n, phase, 1.00, 1.00))
        self._param.dihedral_types[self._dih_type] = dih_tlist
