import simtk.openmm as mm
from torsionfit import TorsionScanSet as ScanSet
import torsionfit.TorsionFitModel as Model
from torsionfit import sqlite_plus
from pymc import MCMC
from parmed.charmm import CharmmParameterSet

param = CharmmParameterSet('../../../../../param/top_all36_cgenff.rtf',
                           '../../../../../../data/charmm_ff/par_all36_cgenff.prm')
structure = '../../../../../structure/butane_charge_off.psf'
scan = '../../../../../torsion_scans/DFT_b3lyp/butane_scan_b3lyp_360.log'

param.dihedral_types[('CG331', 'CG321', 'CG321', 'CG331')][1].phi_k=0
param.dihedral_types[('CG331', 'CG321', 'CG321', 'CG331')][0].phi_k=0
param.dihedral_types[('HGA3', 'CG331', 'CG321', 'HGA2')][0].phi_k=0
param.dihedral_types[('HGA2', 'CG321', 'CG331', 'HGA3')][0].phi_k=0
param.dihedral_types[('HGA3', 'CG331', 'CG321', 'CG321')][0].phi_k=0
param.dihedral_types[('CG321', 'CG321', 'CG331', 'HGA3')][0].phi_k=0
param.dihedral_types[('HGA2', 'CG321', 'CG321', 'HGA2')][0].phi_k=0
param.dihedral_types[('CG331', 'CG321', 'CG321', 'HGA2')][0].phi_k=0

butane_scan = ScanSet.parse_psi4(scan, structure)

platform = mm.Platform.getPlatformByName('Reference')

model = Model.TorsionFitModelEliminatePhase(param, butane_scan, platform=platform, decouple_n=True
                                            param_to_opt=[('CG331', 'CG321', 'CG321', 'CG331'),
                                                          ('HGA3', 'CG331', 'CG321', 'HGA2'),
                                                          ('HGA3', 'CG331', 'CG321', 'CG321'),
                                                          ('HGA2', 'CG321', 'CG321', 'HGA2'),
                                                          ('CG331', 'CG321', 'CG321', 'HGA2')])

sampler = MCMC(model.pymc_parameters, db=sqlite_plus, dbname='butane_360_all_mult_on_charge_off_init_zero.db', verbose=5)
sampler.sample(100000)