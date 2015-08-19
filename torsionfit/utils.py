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
