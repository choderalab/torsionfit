""" Test fragmentation """

__author__ = 'Chaya D. Stern'

from torsionfit.tests.utils import get_fun, has_openeye
from torsionfit.qmscan import fragment
from openmoltools.openeye import get_charges, smiles_to_oemol
import unittest
if has_openeye:
    import openeye.oechem as oechem


mol = smiles_to_oemol('CN(C)C/C=C/C(=O)NC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C=C3)F)Cl)O[C@H]4CCOC4')
charged = get_charges(mol, keep_confs=1)


class TestFragments(unittest.TestCase):

    @unittest.skipUnless(has_openeye, "Cannot test without OpenEye")
    def test_tag_funcgroup(self):
        """ Test tag functional groups """
        tagged_funcgroups = fragment._tag_fgroups(charged)
        self.assertEquals(len(tagged_funcgroups), 3)
        atom_idx = tagged_funcgroups['amide_0'][0].pop()
        atom = charged.GetAtom(oechem.OEHasAtomIdx(atom_idx))
        fgroup = atom.GetData('fgroup')
        self.assertEquals('amide_0', fgroup)

    @unittest.skipUnless(has_openeye, "Cannot test without OpenEye")
    def test_tag_rings(self):
        """ Test tag rings"""
        tagged_rings = fragment._tag_rings(charged)
        self.assertEquals(len(tagged_rings), 3)
        atom_idx = tagged_rings[1][0].pop()
        atom = charged.GetAtom(oechem.OEHasAtomIdx(atom_idx))
        ringsystem = atom.GetData('ringsystem')
        self.assertEquals(1, ringsystem)

    @unittest.skipUnless(has_openeye, "Cannot test without OpenEye")
    def test_union(self):
        """ Test combining ring and functional group"""
        tagged_funcgroup = fragment._tag_fgroups(charged)
        tagged_rings = fragment._tag_rings(charged)
        tagged_funcgroup2 = fragment._ring_fgroup_union(charged, tagged_fgroups=tagged_funcgroup,
                                                        tagged_rings=tagged_rings)
        self.assertNotEqual(len(tagged_funcgroup['ether_1'][0]), len(tagged_funcgroup2['ether_1'][0]))

    @unittest.skipUnless(has_openeye, "Cannot test without OpenEye")
    def test_tag_molecule(self):
        """Test tag molecule"""
        tagged_rings, tagged_funcgroup = fragment.tag_molecule(charged)
        atom_idx = tagged_rings[1][0].pop()
        atom = charged.GetAtom(oechem.OEHasAtomIdx(atom_idx))
        ringsystem = atom.GetData('ringsystem')
        self.assertEquals(1, ringsystem)
        atom_idx = tagged_funcgroup['amide_0'][0].pop()
        atom = charged.GetAtom(oechem.OEHasAtomIdx(atom_idx))
        funcgroup = atom.GetData('fgroup')
        self.assertEquals('amide_0', funcgroup)

    @unittest.skipUnless(has_openeye, "Cannot test without OpenEye")
    def test_is_ortho(self):
        """Test if substituent is ortho to rotatable bond"""
        rot_bond = charged.GetBond(oechem.OEHasBondIdx(8))
        next_bond = charged.GetBond(oechem.OEHasBondIdx(10))
        ortho_bond = charged.GetBond(oechem.OEHasBondIdx(30))
        self.assertTrue(fragment._is_ortho(bond=ortho_bond, rot_bond=rot_bond, next_bond=next_bond))
        other_bond = charged.GetBond(oechem.OEHasBondIdx(31))
        self.assertFalse(fragment._is_ortho(bond=other_bond, rot_bond=rot_bond, next_bond=next_bond))







