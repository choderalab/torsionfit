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
