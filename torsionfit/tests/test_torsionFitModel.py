"""Test TorsionFitModel"""

from torsionfit.tests.utils import get_fun
import torsionfit.database.qmdatabase as qmdb
from torsionfit.model import TorsionFitModel
import torsionfit.parameters as par
from pymc import MCMC
import glob
import pymc
from numpy.testing import assert_
from cclib.parser import Gaussian
from cclib.parser.utils import convertor
from numpy.testing import TestCase
from parmed.charmm import CharmmParameterSet
import unittest


try:
    from simtk.openmm import app
    import simtk.openmm as mm
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
structure = get_fun('butane.psf')
logfiles = get_fun('MP2_torsion_scan/')
frag = qmdb.parse_psi4_out(logfiles, structure)
frag = frag.remove_nonoptimized()
to_optimize = [('CG331', 'CG321', 'CG321', 'CG331')]
model = TorsionFitModel(param=param, frags=frag, param_to_opt=to_optimize)


class TestFitModel(unittest.TestCase):
    """ Tests pymc model"""

    def test_pymc_model(self):
        """ Tests sampler """

        sampler = MCMC(model.pymc_parameters)
        self.assert_(isinstance(model, TorsionFitModel))
        self.assert_(isinstance(sampler, pymc.MCMC))

        sampler.sample(iter=1)

    def test_update_param_continuous(self):
        """ Tests that update parameter updates the reverse dihedral too in continuous  """

        param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        structure = get_fun('butane.psf')
        logfiles = get_fun('MP2_torsion_scan/')
        frag = qmdb.parse_psi4_out(logfiles, structure)
        frag = frag.remove_nonoptimized()
        to_optimize = [('CG331', 'CG321', 'CG321', 'CG331')]

        model = TorsionFitModel(param=param, frags=frag, param_to_opt=to_optimize, sample_phase=True,
                                continuous_phase=True)
        torsion = model.parameters_to_optimize[0]
        torsion_reverse = tuple(reversed(torsion))
        self.assertEqual(param.dihedral_types[torsion], param.dihedral_types[torsion_reverse])

    def test_update_param(self):
        """ Tests that update parameter updates the reverse dihedral too """

        model.update_param(param)
        torsion = model.parameters_to_optimize[0]
        torsion_reverse = tuple(reversed(torsion))
        self.assertEqual(param.dihedral_types[torsion], param.dihedral_types[torsion_reverse])

    def test_update_param_struct(self):
        """ Tests that update parameter updates assigned parameters in the structure """

        model.update_param(param)
        torsion = frag.structure.dihedrals[0]
        self.assertEqual(torsion.type, param.dihedral_types[(torsion.atom1.type, torsion.atom2.type,
                                                               torsion.atom3.type, torsion.atom4.type)])

    def test_update_param_struct_cont(self):
        """ Tests that update parameter updates assigned parameters in the structure """
        param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        structure = get_fun('butane.psf')
        logfiles = get_fun('MP2_torsion_scan/')
        frag = qmdb.parse_psi4_out(logfiles, structure)
        frag = frag.remove_nonoptimized()
        to_optimize = [('CG331', 'CG321', 'CG321', 'CG331')]

        model = TorsionFitModel(param=param, frags=frag, param_to_opt=to_optimize, sample_phase=True,
                                continuous_phase=True)
        model.update_param(param)
        torsion = frag.structure.dihedrals[0]
        self.assertEqual(torsion.type, param.dihedral_types[(torsion.atom1.type, torsion.atom2.type,
                                                               torsion.atom3.type, torsion.atom4.type)])

    def test_add_missing(self):
        """ Tests that add_missing adds missing terms to parameters_to_optimize """

        par.add_missing(param_list=to_optimize, param=param)
        for i in model.frags[0].structure.dihedrals:
            key = (i.atom1.type, i.atom2.type, i.atom3.type, i.atom4.type)
            key_reverse = tuple(reversed(key))
            if key in model.parameters_to_optimize or key_reverse in model.parameters_to_optimize:
                self.assert_(len(i.type) == 5)

    def test_add_missing_cond(self):
        """ Tests that add_missing adds missing terms to parameters_to_optimize """

        param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        structure = get_fun('butane.psf')
        logfiles = get_fun('MP2_torsion_scan/')
        frag = qmdb.parse_psi4_out(logfiles, structure)
        frag = frag.remove_nonoptimized()
        to_optimize = [('CG331', 'CG321', 'CG321', 'HGA2')]

        model = TorsionFitModel(param=param, frags=frag, param_to_opt=to_optimize, sample_phase=True,
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
        self.assertFalse(model.rj)
        self.assertTrue(('CG331_CG321_CG321_CG331_multiplicity_bitstring' not in model.pymc_parameters))

    def test_rj_on(self):
        """Test reversible jump is on"""

        model = TorsionFitModel(param=param, frags=frag, param_to_opt=to_optimize, rj=True)
        self.assertTrue(model.rj)
        self.assertTrue(('CG331_CG321_CG321_CG331_multiplicity_bitstring' in model.pymc_parameters))




