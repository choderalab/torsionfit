""" Test TorsionScanSet """

from torsionfit.tests.utils import get_fun
import torsionfit.TorsionScanSet as torsionset
from cclib.parser import Gaussian
from cclib.parser.utils import convertor
from numpy.testing import assert_almost_equal, assert_array_equal
from parmed.charmm import CharmmParameterSet, CharmmPsfFile
import unittest

try:
    from simtk.openmm import app
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False


def test_read_logfile():
    """Tests the logfile parser"""
    structure = get_fun('PRL.psf')
    logfiles = [get_fun('PRL.scan2.neg.log'), get_fun('PRL.scan2.pos.log')]
    scan = torsionset.read_scan_logfile(logfiles, structure)
    scan1 = Gaussian(logfiles[0])
    scan2 = Gaussian(logfiles[1])
    data1 = scan1.parse()
    data2 = scan2.parse()
    assert_almost_equal(scan.xyz[:40]*10, data1.atomcoords, decimal=6)
    assert_almost_equal(scan.xyz[40:]*10, data2.atomcoords, decimal=6)


def test_converged_structures():
    """Tests that non converged coordinates in Gaussian log files are discarded"""
    structure = get_fun('MPR.psf')
    scan = torsionset.read_scan_logfile(get_fun('MPR.scan1.pos.log'), structure)
    log = Gaussian(get_fun('MPR.scan1.pos.log'))
    data = log.parse()
    assert_array_equal(data.atomcoords.shape, scan.xyz.shape)
    converted = convertor(data.scfenergies, "eV", "kJmol-1") - \
                min(convertor(data.scfenergies[:len(data.atomcoords)], "eV", "kJmol-1"))
    assert_array_equal(converted[:47], scan.qm_energy)


def test_extract_geom_opt():
    """Tests extraction of optimized geometry"""
    structure = get_fun('MPR.psf')
    scan = torsionset.read_scan_logfile([get_fun('MPR.scan1.pos.log'), get_fun('MPR.scan1.neg.log')], structure)
    scan.extract_geom_opt()


class TestScanSet(unittest.TestCase):
    """Tests ScanSet"""

    def test_to_optimize(self):
        """Tests generate to_optimize list"""
        param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        structure = (get_fun('PRL.psf'))
        stream = (get_fun('PRL.str'))
        model_param_to_optimize = torsionset.to_optimize(param, stream)
        scan_set = torsionset.read_scan_logfile([get_fun('PRL.scan2.neg.log'), get_fun('PRL.scan2.pos.log')], structure)
        scan_set.get_params(model_param_to_optimize, param)
        self.assertItemsEqual(model_param_to_optimize, scan_set.to_optimize)

    def test_multiple_to_optimize(self):
        """ Tests to optimize lists for multiple fragments"""
        param1 = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        param2 = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
        structs = [get_fun('PRL.psf'), get_fun('MPR.psf')]
        stream = [get_fun('PRL.str'), get_fun('MPR.str')]
        model_param_to_optimize = torsionset.to_optimize(param1, stream)
        mpr_scan = torsionset.read_scan_logfile(get_fun('MPR.scan1.pos.log'), structs[1])
        mpr_scan.get_params(model_param_to_optimize, param1)
        mpr_to_optimize = torsionset.to_optimize(param2, stream[-1])
        self.assertItemsEqual(mpr_to_optimize, mpr_scan.to_optimize)