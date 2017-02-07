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

