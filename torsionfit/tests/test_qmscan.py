""" Tests generating files for qm torsion scan """

import unittest
from torsionfit.tests.utils import get_fn, has_openeye
import torsionfit.qmscan.torsion_scan as qmscan
import tempfile
import os
from fnmatch import fnmatch
import shutil


class TestQmscan(unittest.TestCase):

    @unittest.skipUnless(has_openeye, 'Cannot test without openeye')
    def test_generat_torsions(self):
        """ Tests finding torsion to drive """
        mol = get_fn('butane.pdb')
        outfile_path = tempfile.mkdtemp()[1]
        qmscan.generate_torsions(mol=mol, path=outfile_path, interval=30)
        input_files = []
        pattern = '*.pdb'
        for path, subdir, files in os.walk(outfile_path):
            for name in files:
                if fnmatch(name, pattern):
                    input_files.append(os.path.join(path, name))

        contents = open(input_files[0]).read()
        pdb = get_fn('butane_10_7_4_3_0.pdb')
        compare_contents = open(pdb).read()
        self.assertEqual(contents, compare_contents )

        shutil.rmtree(outfile_path)

    def test_generate_input(self):
        """Test generate psi4 input files"""
        root = get_fn('torsion_scan/10_7_4_3')
        qmscan.generate_scan_input(root, 'pdb', 'butane', ['MP2'], ['aug-cc-pvtz'], symmetry='C1')

        contents = open(get_fn('torsion_scan/10_7_4_3/0/butane_10_7_4_3_0.dat')).read()
        compare_content = open(get_fn('butane_10_7_4_3_0.dat')).read()
        self.assertEqual(contents, compare_content)


#
