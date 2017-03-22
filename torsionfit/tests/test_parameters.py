""" Test parameters """

__author__ = 'Chaya D. Stern'

from torsionfit.tests.utils import get_fun
from torsionfit.backends import sqlite_plus
import torsionfit.parameters as par
from parmed.charmm.parameters import CharmmParameterSet
import unittest


class TestParameters(unittest.TestCase):
    """ Test parameters.py """

    def test_add_missing(self):
        """ Test add missing parameters """
        param = CharmmParameterSet(get_fun('par_all36_cgenff.prm'), get_fun('top_all36_cgenff.rtf'))
        torsion = ('CG331', 'CG321', 'CG321', 'CG331')
        self.assertEqual(len(param.dihedral_types[torsion]), 2)
        par.add_missing(param=param, param_list=torsion)
        self.assertEqual(len(param.dihedral_types[torsion]), 5)
        par.add_missing(param=param, param_list=torsion, sample_n5=True)
        self.assertEqual(len(param.dihedral_types[torsion]), 6)

    def test_set_phase_0(self):
        """ Test set phases to 0 """
        param = CharmmParameterSet(get_fun('par_all36_cgenff.prm'), get_fun('top_all36_cgenff.rtf'))
        torsion = ('CG331', 'CG321', 'CG321', 'CG331')
        self.assertEqual(param.dihedral_types[torsion][1].phase, 180.0)
        par.set_phase_0(torsion, param)
        self.assertEqual(param.dihedral_types[torsion][1].phase, 0.0)

    def test_param_from_db(self):
        """ Tests parameterizing from database """
        param = CharmmParameterSet(get_fun('par_all36_cgenff.prm'), get_fun('top_all36_cgenff.rtf'))
        torsion = ('CG331', 'CG321', 'CG321', 'CG331')
        db = sqlite_plus.load(get_fun('butane.db'))
        par.add_missing(torsion, param, sample_n5=True)
        self.assertEqual(param.dihedral_types[torsion][0].phi_k, 0.03819)
        self.assertEqual(param.dihedral_types[torsion][1].phi_k, 0.03178)
        self.assertEqual(param.dihedral_types[torsion][2].phi_k, 0.0)
        self.assertEqual(param.dihedral_types[torsion][3].phi_k, 0.0)
        self.assertEqual(param.dihedral_types[torsion][4].phi_k, 0.0)
        self.assertEqual(param.dihedral_types[torsion][5].phi_k, 0.0)

        par.update_param_from_sample(param_list=torsion, param=param, db=db, rj=False)
        self.assertEqual(param.dihedral_types[torsion][0].phi_k, 0.086424)
        self.assertEqual(param.dihedral_types[torsion][1].phi_k, 0.019074)
        self.assertEqual(param.dihedral_types[torsion][2].phi_k, 0.0)
        self.assertEqual(param.dihedral_types[torsion][3].phi_k, -1.834546)
        self.assertEqual(param.dihedral_types[torsion][4].phi_k, -1.86807)
        self.assertEqual(param.dihedral_types[torsion][5].phi_k, -2.860622)





