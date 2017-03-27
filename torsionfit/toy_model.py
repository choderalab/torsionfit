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
import torsionfit.database.qmdatabase as ScanSet
import torsionfit.model as model
from torsionfit.tests.utils import get_fun
import copy


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

    def __init__(self, true_value=None, initial_value=None, n_increments=18, rj=True, sample_phase=False,
                 continuous=False):
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
        if true_value is not None:
            if type(true_value) == DihedralTypeList:
                dih_tlist = true_value
            elif type(true_value) == DihedralType:
                dih_tlist = DihedralTypeList()
                dih_tlist.append(true_value)
        else:
            dih_tlist = self._randomize_dih_param(return_dih=True)
        self.true_value = copy.deepcopy(dih_tlist)
        self._param.dihedral_types[self._dih_type] = dih_tlist

        # parametrize toy
        self._struct.load_parameters(self._param, copy_parameters=False)
        self._struct.positions = self._pdb.positions

        # generate synthetic torsion scan
        self._torsion_scan(n_increments=n_increments)

        # initialize parameter
        if initial_value is not None:
            if type(initial_value) == DihedralTypeList:
                dih_tlist = initial_value
            if type(initial_value) == DihedralType:
                dih_tlist = DihedralTypeList()
                dih_tlist.append(initial_value)
            elif initial_value == 'cgenff':
                dih_tlist = original_torsion
        else:
            dih_tlist = self._randomize_dih_param(return_dih=True)

        self.initial_value = copy.deepcopy(dih_tlist)
        self._param.dihedral_types[self._dih_type] = dih_tlist

        # create torsionfit.TorsionScanSet
        torsions = np.zeros((len(self._positions), 4))
        torsions[:] = [1, 2, 3, 4]
        direction = None
        steps = None
        self.scan_set = ScanSet.QMDataBase(positions=self._positions.value_in_unit(units.nanometers),
                                           topology=self._topology, structure=self._struct, torsions=torsions,
                                           steps=steps, directions=direction,
                                           qm_energies=self.synthetic_energy.value_in_unit(units.kilojoules_per_mole))

        self.model = model.TorsionFitModel(param=self._param, frags=self.scan_set, platform=self._platform,
                                           param_to_opt=[self._dih_type], rj=rj, continuous_phase=continuous,
                                           sample_phase=sample_phase)

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
        phis = np.arange(-np.pi, +np.pi, np.pi/n_increments)
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
            energy = state.getPotentialEnergy()
            self.synthetic_energy[i] = energy + units.Quantity(value=np.random.normal(0, 1.0), unit=units.kilojoules_per_mole)

    def _randomize_dih_param(self, return_dih=False):
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
        if return_dih:
            _dih_tlist = copy.deepcopy(dih_tlist)
            return _dih_tlist

    def save_torsion_values(self, filename=None):
        """

        Parameters
        ----------
        filename : str
            file to save torsion values. Default is None. If None, return np.array

        Returns
        -------
        dih_param: np.array
            array containing torison parameters
        """
        dih_param = np.ones(shape=(2, 6, 3)) * np.nan
        for i, dih in enumerate(self.true_value):
            dih_param[0, i, 0] = dih.per
            dih_param[0, i, 1] = dih.phi_k
            dih_param[0, i, 2] = dih.phase

        for j, dih in enumerate(self.initial_value):
            dih_param[1, j, 0] = dih.per
            dih_param[1, j, 1] = dih.phi_k
            dih_param[1, j, 2] = dih.phase

        if filename is not None:
            np.save(filename, dih_param)
        else:
            return dih_param

    @staticmethod
    def from_dih_params(filename=None, dih_params=None, rj=False, continuous=False, n_increments=13, sample_phase=False):
        """

        Parameters
        ----------
        filename : str
            name of file with serialized dihedral parameters
        rj : bool
            Flag if using reversible jump. Default it True
        continuous : bool
            Flag if sampling continuous phase. Default is False
        n_increments : int
            incermentation of torsion drive
        sample_phase : bool
            Flag if sampling phase. Default is False (K is allowed to go negative when sample_phase is False)

        Returns
        -------
        ToyModel with true and initial value from saved file.

        """
        if filename is None and dih_params is None:
            msg = 'You must provide either an npy file or a numpy array with true and initial values for the toy model'
            raise Exception(msg)
        if filename is not None:
            dih_params = np.load(filename)
        dih_tlist_true = DihedralTypeList()
        dih_tlist_init = DihedralTypeList()

        true = dih_params[0]
        init = dih_params[1]

        for dih in true:
            if not np.isnan(dih[0]):
                dih_tlist_true.append(DihedralType(per=dih[0], phi_k=dih[1], phase=dih[2]))

        for dih in init:
            if not np.isnan(dih[0]):
                dih_tlist_init.append(DihedralType(per=dih[0], phi_k=dih[1], phase=dih[2]))

        return ToyModel(true_value=dih_tlist_true, initial_value=dih_tlist_init, rj=rj, continuous=continuous,
                        n_increments=n_increments, sample_phase=sample_phase)



