import numpy as np

def RMSE(scanSet, db):
    '''

    :param model: TorsionScanSet
    :param db: pymc database
    :return: numpy array of rmse
    '''

    N = len(scanSet.qm_energy)
    errors = np.zeros(len(db.trace('mm_energy')[:]))
    for i, energy in enumerate(db.trace('mm_energy')[:]):
        rmse = np.linalg.norm(energy - scanSet.qm_energy)/np.sqrt(N)
        errors[i] = rmse
    return errors


def get_multiplicity_trace(key, db):
    """
    This function takes a bitstring key and data base and returns a dictionary of [0,1] traces for each multiplicity
    term. 0 when that multiplicity is off and 1 when it's on
    :param key: str (A_B_C_D_multiplicity_bitstring)
    :param db: pymc sqlite_plus database
    :return: dictionary of multiplicity mapped to lists of 0 and 1

    """
    multiplicities = (1, 2, 3, 4, 6)
    multiplicity_trace = {}
    for m in multiplicities:
        multiplicity_trace[str(m)] = []
        for i in db.trace(key)[:]:
            if 2**(m-1) & int(i):
                multiplicity_trace[str(m)].append(1)
            else:
                multiplicity_trace[str(m)].append(0)
    return multiplicity_trace


def get_sampled_torsions(db):
    """
    This function returns a list of torsions that were sampled in the database in the form of A_B_C_D
    :param db: pymc sqlit_plus database

    Returns:
    list of torsions that were sampled

    """
    torsions = []
    for key in db.getstate()['stochastics']:
        key_split = key.split('_')
        if key_split[-1] == 'bitstring':
            name = key_split[0] + '_' + key_split[1] + '_' + key_split[2] + '_' + key_split[3]
            torsions.append(name)
    return torsions


def param_from_db(param, db, i=-1, decouple_n=False, phase=False):
    """
    This function parameterizes sampled torsion with values of sample i in database
     param: parmed.charmm.parameterset
     db: sqlit_plus database
     i: int, sample to use
     decouple_n: flag if multiplicities were sampled. (When true, the n was decoupled and it wasn't sampled)
         Default False
     phase: flag is phases were sampled.

    Returns:

    """
    torsions = get_sampled_torsions(db)
    for torsion_name in torsions:
        multiplicity_bitstring = int(db.trace(torsion_name + '_multiplicity_bitstring')[i])
        t = tuple(torsion_name.split('_'))
        reverse_t = tuple(reversed(t))
        for i in range(len(param.dihedral_types[t])):
            m = int(param.dihedral_types[t][i].per)
            multiplicity_bitmask = 2 ** (m - 1)  # multiplicity bitmask
            if (multiplicity_bitstring & multiplicity_bitmask) or decouple_n:
                if m == 5:
                    continue
                k = torsion_name + '_' + str(m) + '_K'
                sample = db.trace(k)[i]
                param.dihedral_types[t][i].phi_k = sample
                param.dihedral_types[reverse_t][i].phi_k = sample
                if phase:
                    p = torsion_name + '_' + str(m) + '_Phase'
                    sample = db.trace(p)[i]
                    param.dihedral_types[t][i].phase = sample
                    param.dihedral_types[reverse_t][i].phase = sample
            else:
                # This torsion periodicity is disabled.
                param.dihedral_types[t][i].phi_k = 0
                param.dihedral_types[reverse_t][i].phi_k = 0