""" Test TorsionScanSet """

from torsionfit.testing.testing import get_fun
import torsionfit.TorsionScanSet as torsionset
from cclib.parser import Gaussian
from numpy.testing import assert_almost_equal

try:
    from simtk.openmm import app
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False


def test_read_logfile():
    structure = get_fun('PRL.psf')
    logfiles = [get_fun('PRL.scan2.neg.log'), get_fun('PRL.scan2.pos.log')]
    scan = torsionset.read_scan_logfile(logfiles, structure)
    scan1 = Gaussian(logfiles[0])
    scan2 = Gaussian(logfiles[1])
    data1 = scan1.parse()
    data2 = scan2.parse()
    assert_almost_equal(scan.xyz[:40]*10, data1.atomcoords, decimal=6)
    assert_almost_equal(scan.xyz[40:]*10, data2.atomcoords, decimal=6)

