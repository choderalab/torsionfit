__author__ = 'Chaya D. Stern'

import numpy as np
import logging
import sys

verbose = False

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


def logger(name='torsionFit', pattern='%(asctime)s %(levelname)s %(name)s: %(message)s',
           date_format='%H:%M:%S', handler=logging.StreamHandler(sys.stdout)):
    """
    Retrieves the logger instance associated to the given name
    :param name: The name of the logger instance
    :param pattern: The associated pattern
    :param date_format: The date format to be used in the pattern
    :param handler: The logging handler
    :return: The logger
    """
    _logger = logging.getLogger(name)
    _logger.setLevel(log_level(verbose))

    if not _logger.handlers:
        formatter = logging.Formatter(pattern, date_format)
        handler.setFormatter(formatter)
        handler.setLevel(log_level(verbose))
        _logger.addHandler(handler)
        _logger.propagate = False
    return _logger


def log_level(verbose=verbose):
    if verbose:
        return logging.DEBUG
    else:
        return logging.INFO


def param_from_db(param, db, i=-1, decouple_n=False, phase=False, n_5=True):
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
