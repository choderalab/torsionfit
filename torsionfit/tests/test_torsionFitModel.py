"""Test TorsionFitModel"""

from torsionfit.tests.utils import get_fn
import torsionfit.database.qmdatabase as qmdb
from torsionfit.model_omm import TorsionFitModel as TorsionFitModelOMM
from torsionfit.model import TorsionFitModel
from torsionfit.backends import sqlite_plus

import torsionfit.parameters as par
from pymc import MCMC
import pymc
from parmed.charmm import CharmmParameterSet
import unittest

try:
    from simtk.openmm import app
    import simtk.openmm as mm
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

param = CharmmParameterSet(get_fn('top_all36_cgenff.rtf'), get_fn('par_all36_cgenff.prm'))
structure = get_fn('butane.psf')
logfiles = get_fn('MP2_torsion_scan/')
frag = qmdb.parse_psi4_out(logfiles, structure)
frag = frag.remove_nonoptimized()
to_optimize = [('CG331', 'CG321', 'CG321', 'CG331')]
model_omm = TorsionFitModelOMM(param=param, frags=frag, param_to_opt=to_optimize, init_random=False)


class TestFitModel(unittest.TestCase):
    """ Tests pymc model"""

    def test_pymc_model(self):
        """ Tests sampler """

        sampler = MCMC(model_omm.pymc_parameters)
        self.assert_(isinstance(model_omm, TorsionFitModelOMM))
        self.assert_(isinstance(sampler, pymc.MCMC))

        sampler.sample(iter=1)

    def test_update_param_continuous(self):
        """ Tests that update parameter updates the reverse dihedral too in continuous  """

        param = CharmmParameterSet(get_fn('top_all36_cgenff.rtf'), get_fn('par_all36_cgenff.prm'))
        structure = get_fn('butane.psf')
        logfiles = get_fn('MP2_torsion_scan/')
        frag = qmdb.parse_psi4_out(logfiles, structure)
        frag = frag.remove_nonoptimized()
        to_optimize = [('CG331', 'CG321', 'CG321', 'CG331')]

        model = TorsionFitModelOMM(param=param, frags=frag, param_to_opt=to_optimize, sample_phase=True,
                                   continuous_phase=True)
        torsion = model.parameters_to_optimize[0]
        torsion_reverse = tuple(reversed(torsion))
        self.assertEqual(param.dihedral_types[torsion], param.dihedral_types[torsion_reverse])

    def test_update_param(self):
        """ Tests that update parameter updates the reverse dihedral too """

        par.update_param_from_sample(model_omm.parameters_to_optimize, param, model=model_omm, rj=model_omm.rj,
                                     phase=model_omm.sample_phase, n_5=model_omm.sample_n5, model_type='openmm')
        torsion = model_omm.parameters_to_optimize[0]
        torsion_reverse = tuple(reversed(torsion))
        self.assertEqual(param.dihedral_types[torsion], param.dihedral_types[torsion_reverse])

    def test_update_param_struct(self):
        """ Tests that update parameter updates assigned parameters in the structure """

        par.update_param_from_sample(model_omm.parameters_to_optimize, param, model=model_omm, rj=model_omm.rj,
                                     phase=model_omm.sample_phase, n_5=model_omm.sample_n5, model_type='openmm')
        torsion = frag.structure.dihedrals[0]
        self.assertEqual(torsion.type, param.dihedral_types[(torsion.atom1.type, torsion.atom2.type,
                                                               torsion.atom3.type, torsion.atom4.type)])

    def test_update_param_struct_cont(self):
        """ Tests that update parameter updates assigned parameters in the structure """
        param = CharmmParameterSet(get_fn('top_all36_cgenff.rtf'), get_fn('par_all36_cgenff.prm'))
        structure = get_fn('butane.psf')
        logfiles = get_fn('MP2_torsion_scan/')
        frag = qmdb.parse_psi4_out(logfiles, structure)
        frag = frag.remove_nonoptimized()
        to_optimize = [('CG331', 'CG321', 'CG321', 'CG331')]

        model = TorsionFitModelOMM(param=param, frags=frag, param_to_opt=to_optimize, sample_phase=True,
                                continuous_phase=True)
        par.update_param_from_sample(model.parameters_to_optimize, param, model=model, rj=model.rj,
                                     phase=model.sample_phase, n_5=model.sample_n5, model_type='openmm')
        torsion = frag.structure.dihedrals[0]
        self.assertEqual(torsion.type, param.dihedral_types[(torsion.atom1.type, torsion.atom2.type,
                                                               torsion.atom3.type, torsion.atom4.type)])

    def test_add_missing(self):
        """ Tests that add_missing adds missing terms to parameters_to_optimize """

        par.add_missing(param_list=to_optimize, param=param)
        for i in model_omm.frags[0].structure.dihedrals:
            key = (i.atom1.type, i.atom2.type, i.atom3.type, i.atom4.type)
            key_reverse = tuple(reversed(key))
            if key in model_omm.parameters_to_optimize or key_reverse in model_omm.parameters_to_optimize:
                self.assert_(len(i.type) == 5)

    def test_add_missing_cond(self):
        """ Tests that add_missing adds missing terms to parameters_to_optimize """

        param = CharmmParameterSet(get_fn('top_all36_cgenff.rtf'), get_fn('par_all36_cgenff.prm'))
        structure = get_fn('butane.psf')
        logfiles = get_fn('MP2_torsion_scan/')
        frag = qmdb.parse_psi4_out(logfiles, structure)
        frag = frag.remove_nonoptimized()
        to_optimize = [('CG331', 'CG321', 'CG321', 'HGA2')]

        model = TorsionFitModelOMM(param=param, frags=frag, param_to_opt=to_optimize, sample_phase=True,
                                continuous_phase=True)

        par.add_missing(param_list=to_optimize, param=param)
        for i in model.frags[0].structure.dihedrals:
            key = (i.atom1.type, i.atom2.type, i.atom3.type, i.atom4.type)
            key_reverse = tuple(reversed(key))
            if key in model.parameters_to_optimize or key_reverse in model.parameters_to_optimize:
                self.assert_(len(i.type) == 5)

    def test_rj(self):
        """ Test reversible jump is off

        """
        self.assertFalse(model_omm.rj)
        self.assertTrue(('CG331_CG321_CG321_CG331_multiplicity_bitstring' not in model_omm.pymc_parameters))

    def test_rj_on(self):
        """Test reversible jump is on"""

        model = TorsionFitModelOMM(param=param, frags=frag, param_to_opt=to_optimize, rj=True)
        self.assertTrue(model.rj)
        self.assertTrue(('CG331_CG321_CG321_CG331_multiplicity_bitstring' in model.pymc_parameters))

    def test_phase(self):
        """ Test phase on and off

        """
        self.assertFalse(model_omm.sample_phase)
        self.assertFalse('CG331_CG321_CG321_CG331_1_Phase' in model_omm.pymc_parameters)

        model_phase = TorsionFitModelOMM(param=param, frags=frag, param_to_opt=to_optimize, sample_phase=True)
        self.assertTrue(model_phase.sample_phase)
        self.assertTrue('CG331_CG321_CG321_CG331_1_Phase' in model_phase.pymc_parameters)

    def test_residual_energy(self):
        """ Tests that total energy is resonable.
        """
        db = sqlite_plus.load(get_fn('butane_np.sqlite'))

        dih_list = [('CG331', 'CG321', 'CG321', 'CG331'),
                    ('HGA2', 'CG321', 'CG321', 'HGA2'),
                    ('CG331', 'CG321', 'CG321', 'HGA2')]

        par.add_missing(param_list=dih_list, param=param, sample_n5=True)
        par.update_param_from_sample(param_list=dih_list, param=param, db=db, n_5=True, rj=False, model_type='numpy')

        frag.compute_energy(param=param)

        self.assertTrue((frag.delta_energy._value > -0.5).all() and (frag.delta_energy._value < 0.5).all())



