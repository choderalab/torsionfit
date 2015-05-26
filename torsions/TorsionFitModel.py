from pymc import Uniform, DiscreteUniform
import TorsionScanSet
import numpy as np


class TorsionFitModel(object):

    def __init__(self, param, stream, frags):

        self.pymc_parameters = dict()
        self.multiplicity_bitstring = dict()

        multiplicities = [1, 2, 3, 4, 6]
        multiplicity_bitstrings = dict()
        self.parameters_to_optimize = TorsionScanSet.to_optimize(param, stream)
        for p in self.parameters_to_optimize:
            torsion_name = p[0]+'_'+p[1]+'_'+p[2]+'_'+p[3]
            if torsion_name not in multiplicity_bitstrings.keys():
                multiplicity_bitstrings[torsion_name] = 0

            for m in multiplicities:
                name = p[0]+'_'+p[1]+'_'+p[2]+'_'+p[3]+'_' + str(m) + '_K'
                k = Uniform(name, lower=0, upper=20, value=0)
                for i in range(len(param.dihedral_types[p])):
                    if param.dihedral_types[p][i].per == m:
                        multiplicity_bitstrings[torsion_name] += 2**(m-1)
                        k = Uniform(name, lower=0, upper=20, value=param.dihedral_types[p][i].phi_k)
                        break

                self.pymc_parameters[name] = k

                name = p[0]+'_'+p[1]+'_'+p[2]+'_'+p[3]+'_' + str(m) + '_Phase'
                phase = DiscreteUniform(name, lower=0, upper=1, value=0)
                for i in range(len(param.dihedral_types[p])):
                    if param.dihedral_types[p][i].per == m:
                        if param.dihedral_types[p][i].phase == 0:
                            phase = DiscreteUniform(name, lower=0, upper=1, value=0)
                            break
                        if param.dihedral_types[p][i].phase == 3.141592653589793:
                            phase = DiscreteUniform(name, lower=0, upper=1, value=1)
                            break

                self.pymc_parameters[name] = phase

        for torsion_name in multiplicity_bitstrings.keys():
            name = torsion_name + '_multiplicity_bitstring'
            bitstring = DiscreteUniform(name, lower=0, upper=63, value=multiplicity_bitstrings[torsion_name])
            self.multiplicity_bitstring[name] = bitstring

        self.frags = frags

    def update_param(self, param):
        """
        Update param set based on current pymc model parameters.

        :return:
        """
        multiplicities = [1, 2, 3, 4, 6]
        # multiplicity_bitstrings = dict()
        for p in self.parameters_to_optimize:
            torsion_name = p[0]+'_'+p[1]+'_'+p[2]+'_'+p[3]
            # multiplicity_bitstring = multiplicity_bitstrings[torsion_name]
            multiplicity_bitstring = self.multiplicity_bitstring[torsion_name + '_multiplicity_bitstring'].value
            print multiplicity_bitstring

            for i in range(len(param.dihedral_types[p])):
                m = int(param.dihedral_types[p][i].per)
                multiplicity_bitmask = 2**(m-1) # multiplicity bitmask
                print multiplicity_bitmask
                if multiplicity_bitstring & multiplicity_bitmask:
                    print multiplicity_bitstring
                    print multiplicity_bitmask
                    print "multiplicity = %s" % m
                    print "True"
                    if m == 5:
                        continue
                    k = torsion_name + '_' + str(m) + '_K'
                    phase = torsion_name + '_' + str(m) + '_Phase'
                    # pymc_variable = getattr(self, name)
                    pymc_variable = self.pymc_parameters[k]
                    param.dihedral_types[p][i].phi_k = pymc_variable.value
                    pymc_variable = self.pymc_parameters[phase]
                    param.dihedral_types[p][i].phase = pymc_variable.value
                else:
                    # This torsion periodicity is disabled.
                    param.dihedral_types[p][i].phi_k = 0





