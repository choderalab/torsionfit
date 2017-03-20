"""
Useful functions for manipulating parameters in a Parmed CharmmParameterSet.

"""
__author__ = 'Chaya D. Stern'

from parmed.topologyobjects import DihedralType


def add_missing(param_list, param, sample_n5=False, sample_phase=True):

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


def param_from_db(param, db, i=-1, decouple_n=False, phase=False, n_5=True):
    """
    This function parameterizes sampled torsion with values of sample i in database. The modifications are in place.

    parameters:
    -----------
     param: parmed.charmm.parameterset
     db: sqlit_plus database
     i: int, sample to use
     decouple_n: flag if multiplicities were sampled. (When true, the n was decoupled and it wasn't sampled)
         Default False
     phase: flag if phases were sampled.
     n_5: bool
        Flag if multiplicity of 5 was sampled and also needs to be modified. Default is True.
    """
    torsions = db.get_sampled_torsions()
    for torsion_name in torsions:
        multiplicity_bitstring = int(db.trace(torsion_name + '_multiplicity_bitstring')[i])
        t = tuple(torsion_name.split('_'))
        reverse_t = tuple(reversed(t))
        for n in range(len(param.dihedral_types[t])):
            m = int(param.dihedral_types[t][n].per)
            multiplicity_bitmask = 2 ** (m - 1)  # multiplicity bitmask
            if (multiplicity_bitstring & multiplicity_bitmask) or decouple_n:

                if m == 5 and not n_5:
                    continue
                k = torsion_name + '_' + str(m) + '_K'
                sample = db.trace(k)[i]
                param.dihedral_types[t][n].phi_k = sample
                param.dihedral_types[reverse_t][n].phi_k = sample
                if phase:
                    p = torsion_name + '_' + str(m) + '_Phase'
                    sample = db.trace(p)[i]
                    param.dihedral_types[t][n].phase = sample
                    param.dihedral_types[reverse_t][n].phase = sample
            else:
                # This torsion periodicity is disabled.
                param.dihedral_types[t][n].phi_k = 0
                param.dihedral_types[reverse_t][n].phi_k = 0


# def update_param(self, param):
#     """
#     Update param set based on current pymc model parameters.
#
#     :mol: torsionfit.TorsionScanSet
#
#     :return: updated torsionfit.TorsionScanSet parameters based on current TorsionFitModel parameters
#     """
#
#     for p in self.parameters_to_optimize:
#         torsion_name = p[0] + '_' + p[1] + '_' + p[2] + '_' + p[3]
#         multiplicity_bitstring = self.pymc_parameters[torsion_name + '_multiplicity_bitstring'].value
#         reverse_p = tuple(reversed(p))
#         for i in range(len(param.dihedral_types[p])):
#             m = int(param.dihedral_types[p][i].per)
#             multiplicity_bitmask = 2 ** (m - 1)  # multiplicity bitmask
#             if (multiplicity_bitstring & multiplicity_bitmask) or self.decouple_n:
#                 if m == 5:
#                     continue
#                 k = torsion_name + '_' + str(m) + '_K'
#                 phase = torsion_name + '_' + str(m) + '_Phase'
#                 pymc_variable = self.pymc_parameters[k]
#                 param.dihedral_types[p][i].phi_k = pymc_variable.value
#                 param.dihedral_types[reverse_p][i].phi_k = pymc_variable.value
#                 pymc_variable = self.pymc_parameters[phase].value
#                 if pymc_variable == 1:
#                     param.dihedral_types[p][i].phase = 180
#                     param.dihedral_types[reverse_p][i].phase = 180
#                     break
#
#                 if pymc_variable == 0:
#                     param.dihedral_types[p][i].phase = 0
#                     param.dihedral_types[reverse_p][i].phase = 0
#                     break
#             else:
#                 # This torsion periodicity is disabled.
#                 param.dihedral_types[p][i].phi_k = 0
#                 param.dihedral_types[reverse_p][i].phi_k = 0







