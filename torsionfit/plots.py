"""
Plotting module for torsionfit

"""

__author__ = 'Chaya D. Stern'

import matplotlib.pyplot as plt
import numpy as np


def marg_mult(model, db, burn, filename):
    """
    generates histogram for marginal distribution of posterior multiplicities.

    :param model: TorsionFitModel
    :param db: pymc.database for model
    :param burn: int. number of steps to skip
    :param filename: filename for plot to save
    """
    mult_bitstring = []
    for i in model.pymc_parameters.keys():
        if i.split('_')[-1] == 'bitstring':
            mult_bitstring.append(i)

    bitmask = [1, 2, 3, 4, 6]
    multiplicities = np.zeros((len(mult_bitstring), 1000000, 5))

    for m, torsion in enumerate(mult_bitstring):
        for i, j in enumerate(db.trace('%s' % torsion)[burn:]):
            for k, l in enumerate(bitmask):
                if 2**(l-1) & int(j):
                    multiplicities[m][i][k] = 1

    plt.matshow(multiplicities.sum(1), cmap='jet',  extent=[0, 5, 0, 20]), plt.colorbar()
    plt.yticks([])
    plt.xlabel('multiplicity term')
    plt.ylabel('torsion')
    plt.savefig(filename)

