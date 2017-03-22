"""
Useful functions for manipulating parameters in a Parmed CharmmParameterSet.

"""
__author__ = 'Chaya D. Stern'

from parmed.topologyobjects import DihedralType
from torsionfit.utils import logger


def add_missing(param_list, param, sample_n5=False):

    """
    Update param set with missing multiplicities. The modifications are in place.

    Parameters
    ----------
    param_list : list of tuples
        list of parameters to add missing. Format (A, B, C, D)
    param : parmed.charmm.CharmmParameterSet
        parameter set to add missing parameters.
    sample_n5 : bool
        Flag if multiplicity of 5 should be added. Default False.
    """
    multiplicities = [1, 2, 3, 4, 6]
    if sample_n5:
        multiplicities = [1, 2, 3, 4, 5, 6]
    if type(param_list) is not list:
        param_list = [param_list]
    for p in param_list:
        reverse = tuple(reversed(p))
        per = []
        for i in range(len(param.dihedral_types[p])):
            per.append(param.dihedral_types[p][i].per)
            per.append(param.dihedral_types[reverse][i].per)
        for j in multiplicities:
            if j not in per:
                param.dihedral_types[p].append(DihedralType(0, j, 0))
                param.dihedral_types[reverse].append(DihedralType(0, j, 0))


def set_phase_0(param_list, param):
    """
    Set all phase angles to 0. Modifications are in place.

    Parameters
    ----------
    param_list : list of tuples.
        list of parameters to set phase to 0. Format (A, B, C, D)
    param : parmed CharmmParameterSet

    """
    if type(param_list) is not list:
        param_list = [param_list]
    for p in param_list:
        reverse_p = tuple(reversed(p))
        for i in range(len(param.dihedral_types[p])):
            param.dihedral_types[p][i].phase = 0
            param.dihedral_types[reverse_p][i].phase = 0


def update_param_from_sample(param_list, param, db=None, model=None, i=-1, rj=False, phase=False, n_5=True, continuous=False):
    """
    This function parameterizes sampled torsion with values of sample i in database or current value in pymc model.
    The modifications are in place.

    parameters:
    -----------
     param_list: list
      list of tuples of torsions being sampled [(A, B, C, D), (E, F, G, H)]
     param: parmed.charmm.parameterset
     db: sqlit_plus database or pymc sampler
        default is None
     model: pymc model
        default is None
     i: int, sample to use
        default is -1
     rj: flag if reversible jump is on.
         Default False
     phase: bool
        Flag if phases were sampled. Default is False
     n_5: bool
        Flag if multiplicity of 5 was sampled and also needs to be modified. Default is True.
    """
    logger().debug('updating parameters')
    if type(param_list) is not list:
        param_list = [param_list]
    for t in param_list:
        torsion_name = t[0] + '_' + t[1] + '_' + t[2] + '_' + t[3]
        if rj:
            multiplicity_key = torsion_name + '_multiplicity_bitstring'
            if db is not None:
                multiplicity_bitstring = int(db.trace(multiplicity_key)[i])
            if model is not None:
                multiplicity_bitstring = model.pymc_parameters[multiplicity_key].value
        else:
            multiplicity_bitstring = 65
        reverse_t = tuple(reversed(t))
        for n in range(len(param.dihedral_types[t])):
            m = int(param.dihedral_types[t][n].per)
            logger().debug('Working on {}'.format(m))
            multiplicity_bitmask = 2 ** (m - 1)  # multiplicity bitmask
            if (multiplicity_bitstring & multiplicity_bitmask) or not rj:
                if m == 5 and not n_5:
                    continue
                k = torsion_name + '_' + str(m) + '_K'
                if db is not None:
                    sample = db.trace(k)[i]
                if model is not None:
                    sample = model.pymc_parameters[k].value
                logger().debug('K sample value {}'.format(sample))
                param.dihedral_types[t][n].phi_k = sample
                param.dihedral_types[reverse_t][n].phi_k = sample
                if phase:
                    p = torsion_name + '_' + str(m) + '_Phase'
                    if db is not None:
                        sample = db.trace(p)[i]
                    if model is not None:
                        sample = model.pymc_parameters[p].value
                    if not continuous:
                        logger().debug('Not continuous')
                        if sample == 1:
                            sample = 180.0
                    logger().debug('Phase sample value {}'.format(sample))
                    param.dihedral_types[t][n].phase = sample
                    param.dihedral_types[reverse_t][n].phase = sample
            else:
                # This torsion periodicity is disabled.
                logger().debug('Turning off {}'.format(m))
                param.dihedral_types[t][n].phi_k = 0
                param.dihedral_types[reverse_t][n].phi_k = 0
