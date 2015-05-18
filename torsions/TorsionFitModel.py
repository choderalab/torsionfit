from pymc import Uniform, DiscreteUniform
import TorsionScanSet
import numpy as np


class TorsionFitModel(object):

    def __init__(self, param, stream):
        multiplicities = ['1', '2', '3', '4', '6']
        for p in TorsionScanSet.to_optimize(param, stream):
            for m in multiplicities:
                name = p[0]+'_'+p[1]+'_'+p[2]+'_'+p[3]+'_' + m + '_K'
                for i in range(len(param.dihedral_types[p])):
                    if param.dihedral_types[p][i].per == m:
                        k = Uniform(name, lower=0, upper=3, value=param.dihedral_types[p][i].phi_k)
                    else:
                        k = Uniform(name, lower=0, upper=3, value=np.random.uniform(low=0.0, high=3.0))

                setattr(self, name, k)

                name = p[0]+'_'+p[1]+'_'+p[2]+'_'+p[3]+'_' + m + '_Phase'
                for i in range(len(param.dihedral_types[p])):
                    if param.dihedral_types[p][i].per == m:
                        if param.dihedral_types[p][i].phase == 0:
                            value = 0
                        else:
                            value = 1
                        phase = DiscreteUniform(name, lower=0, upper=1, value=value)
                    else:
                        phase = DiscreteUniform(name, lower=0, upper=1, value=np.random.randint(2))

                setattr(self, name, phase)




