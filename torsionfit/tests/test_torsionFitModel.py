"""Test TorsionFitModel"""

from torsionfit.tests.utils import get_fun, get_files
import torsionfit.TorsionScanSet as torsionset
from torsionfit.TorsionFitModel import TorsionFitModel
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

platform = mm.Platform.getPlatformByName('Reference')
param = CharmmParameterSet(get_fun('top_all36_cgenff.rtf'), get_fun('par_all36_cgenff.prm'))
stream = get_fun('PRL.str')
structure = get_fun('PRL.psf')
logfiles = get_files(glob.glob('PRL.scan*'))
frag = torsionset.read_scan_logfile(logfiles, structure)
frag = frag.extract_geom_opt()
model = TorsionFitModel(param, stream, frag, platform=platform)
sampler =MCMC(model.pymc_parameters)

class TestFitModel(unittest.TestCase):
    """ Tests pymc model"""

    def test_pymc_model(self):
        """ Tests sampler """
        self.assert_(isinstance(model, TorsionFitModel))
        self.assert_(isinstance(sampler, pymc.MCMC))
        sampler.sample(iter=1)

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

    def test_add_missing(self):
        """ Tests that add_missing adds missing terms to parameters_to_optimize """

        model.add_missing(param)
        for i in model.frags[0].structure.dihedrals:
            key = (i.atom1.type, i.atom2.type, i.atom3.type, i.atom4.type)
            key_reverse = tuple(reversed(key))
            if key in model.parameters_to_optimize or key_reverse in model.parameters_to_optimize:
                self.assert_(len(i.type) == 5)





