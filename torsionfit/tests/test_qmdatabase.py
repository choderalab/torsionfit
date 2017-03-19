""" Test TorsionScanSet """

from torsionfit.tests.utils import get_fun
import torsionfit.database.qmdatabase as qmdb
from cclib.parser import Gaussian
from cclib.parser.utils import convertor
from numpy.testing import assert_equal, assert_almost_equal, assert_array_equal
from parmed.charmm import CharmmParameterSet
import unittest
import numpy as np


try:
    from simtk.openmm import app
    import simtk.openmm as mm
    import simtk.unit as u
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

# create scan set


def test_read_logfile():
    """Tests the logfile parser"""
    structure = get_fun('PRL.psf')
    logfiles = [get_fun('PRL.scan2.neg.log'), get_fun('PRL.scan2.pos.log')]
    scan = qmdb.parse_gauss(logfiles, structure)
    scan1 = Gaussian(logfiles[0])
    scan2 = Gaussian(logfiles[1])
    data1 = scan1.parse()
    data2 = scan2.parse()
    assert_almost_equal(scan.xyz[:40]*10, data1.atomcoords, decimal=6)
    assert_almost_equal(scan.xyz[40:]*10, data2.atomcoords, decimal=6)


def test_converged_structures():
    """Tests that non converged coordinates in Gaussian log files are discarded"""
    structure = get_fun('MPR.psf')
    scan = qmdb.parse_gauss(get_fun('MPR.scan1.pos.log'), structure)
    log = Gaussian(get_fun('MPR.scan1.pos.log'))
    data = log.parse()
    assert_array_equal(data.atomcoords.shape, scan.xyz.shape)
    converted = convertor(data.scfenergies, "eV", "kJmol-1") - \
                min(convertor(data.scfenergies[:len(data.atomcoords)], "eV", "kJmol-1"))
    assert_array_equal(converted[:47], scan.qm_energy)


def test_extract_geom_opt():
    """Tests extraction of optimized geometry"""
    structure = get_fun('MPR.psf')
    scan = qmdb.parse_gauss([get_fun('MPR.scan1.pos.log'), get_fun('MPR.scan1.neg.log')], structure)
    scan.extract_geom_opt()


class TestScanSet(unittest.TestCase):
    """Tests ScanSet"""

    def test_to_optimize(self):
        """Tests generate to_optimize list"""
        param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        stream = (get_fun('PRL.str'))
        to_optimize = qmdb.to_optimize(param, stream)
        self.assert_(len(to_optimize) == 19)
        reverse = tuple(reversed(to_optimize[0]))
        self.assert_(reverse not in to_optimize)

    def test_parse_psi4_out(self):
        """ Tests psi4 outfile parser"""
        structure = get_fun('butane.psf')
        scan = get_fun('MP2_torsion_scan/')
        butane_scan = qmdb.parse_psi4_out(scan, structure)
        butane_scan = butane_scan.remove_nonoptimized()
        self.assertEqual(butane_scan.n_atoms, 14)
        self.assertEqual(butane_scan.n_chains, 1)
        self.assertEqual(butane_scan.n_residues, 1)
        self.assertEqual(butane_scan.n_frames, 13)
        self.assertEqual(butane_scan.qm_energy.shape, (13,))
        angles = np.arange(0, 370, 30)
        np.testing.assert_equal(butane_scan.angles, angles)
        torsion = np.array([3, 6, 9, 13])
        np.testing.assert_equal(butane_scan.torsion_index[0], torsion)

    def test_remove_nonoptimized(self):
        """ Test remove non_optimized structures """
        structure = get_fun('butane.psf')
        scan = get_fun('MP2_torsion_scan/')
        test_scan = qmdb.parse_psi4_out(scan, structure)
        self.assertEqual(test_scan.n_frames, 14)
        scan_opt = test_scan.remove_nonoptimized()
        self.assertEqual(scan_opt.n_frames, 13)

    def test_to_dataframe(self):
        """ Tests to dataframe """
        structure = get_fun('butane.psf')
        scan = get_fun('MP2_torsion_scan/')
        test_scan = qmdb.parse_psi4_out(scan, structure)
        test_scan.to_dataframe()
        structure = get_fun('MPR.psf')
        scan = qmdb.parse_gauss([get_fun('MPR.scan1.pos.log'), get_fun('MPR.scan1.neg.log')], structure)
        scan.extract_geom_opt()
        scan.to_dataframe(psi4=False)

    def test_compute_energy(self):
        """ Tests compute mm energy"""
        structure = get_fun('butane.psf')
        scan = get_fun('MP2_torsion_scan/')
        test_scan = qmdb.parse_psi4_out(scan, structure)
        param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        self.assertFalse(test_scan._have_mm_energy)
        scan_opt = test_scan.remove_nonoptimized()
        scan_opt.compute_energy(param)
        self.assertTrue(scan_opt._have_mm_energy)
        mm_energy = np.array([22.66381775,  13.11040092,   3.04552792,   7.83767718,
                              15.21426107,   7.0116804 ,   0.        ,   7.00802623,
                              15.21461956,   7.83696204,   3.04525798,  13.10813678,  22.66375837])
        np.testing.assert_almost_equal(scan_opt.mm_energy._value, mm_energy, 4)

    def test_create_context(self):
        """ Test create context """
        structure = get_fun('butane.psf')
        scan = get_fun('MP2_torsion_scan/')
        test_scan = qmdb.parse_psi4_out(scan, structure)
        param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        test_scan.integrator = mm.VerletIntegrator(0.004*u.picoseconds)
        test_scan.create_context(param)

    def test_copy_torsions(self):
        """ Test copy torsions"""
        structure = get_fun('butane.psf')
        scan = get_fun('MP2_torsion_scan/')
        test_scan = qmdb.parse_psi4_out(scan, structure)
        param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        test_scan.compute_energy(param)
        test_scan.copy_torsions()

    def test_mm_from_param_sample(self):
        """"""
        pass

