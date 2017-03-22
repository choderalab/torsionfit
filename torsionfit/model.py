__author__ = 'Chaya D. Stern'

import pymc
import numpy as np
from simtk.unit import kilojoules_per_mole
import torsionfit.database.qmdatabase as TorsionScan
import torsionfit.parameters as par
import warnings


class TorsionFitModel(object):
    """pymc model

    This model only allows a phase angle of 0 but allows force constants to flip signs. If the sign is negative, the
    phase angle will be 180.

    Attributes:
    ----------
    pymc_parameters: dict() of pymc parameters
    parameters_to_optimize: list of tuples (dihedrals to optimize)
    fags: list of TorsionScanSet for fragments
    platform: OpenMM platform to use for potential energy calculations

    """
    def __init__(self, param, frags, stream=None,  platform=None, param_to_opt=None, rj=False, sample_n5=False,
                 continuous_phase=False, sample_phase=False):
        """

        Parameters
        ----------
        param : Parmed CharmmParameterSet
        frags : list of torsionfit.QMDataBase
        stream : str
            Path to CHARMM stream file. Default None.
        platform : openmm.Platform
            Default None.
        param_to_opt : list of tuples of torsions.
            Default None.
        rj : bool
            If True, will use reversible jump to sample Fourier terms. If False, will sample all Ks. Default False
        sample_n5 : bool
            If True, will also sample n=5. Default False
        eliminate_phase : bool
            If True, will not sample phase. Instead, Ks will be able to take on negative values. Default True. If True,
            make sure continuous_phase is also False.
        continuous_phase : bool
            If True, will allow phases to take on any value between 0-180. If False, phase will be a discrete and only
            sample 0 or 180


        Returns
        -------
        pymc model

        """

        if type(frags) != list:
            frags = [frags]

        self.pymc_parameters = dict()
        self.frags = frags
        self.platform = platform
        self.rj = rj
        self.sample_n5 = sample_n5
        self.continuous_phase = continuous_phase
        self.sample_phase = sample_phase
        if param_to_opt:
            self.parameters_to_optimize = param_to_opt
        else:
            self.parameters_to_optimize = TorsionScan.to_optimize(param, stream)

        # Check that options are reasonable
        if not sample_phase and continuous_phase:
            warnings.warn("You can't eliminate phase but have continuous phase. Changing continuous phase to False")
            self.continuous_phase = False

        # set all phases to 0 if eliminate phase is True
        if not self.sample_phase:
            par.set_phase_0(self.parameters_to_optimize, param)

        multiplicities = [1, 2, 3, 4, 6]
        if self.sample_n5:
            multiplicities = [1, 2, 3, 4, 5, 6]
        multiplicity_bitstrings = dict()

        # offset
        for frag in self.frags:
            name = '%s_offset' % frag.topology._residues[0]
            offset = pymc.Uniform(name, lower=-50, upper=50, value=0)
            self.pymc_parameters[name] = offset

        for p in self.parameters_to_optimize:
            torsion_name = p[0] + '_' + p[1] + '_' + p[2] + '_' + p[3]

            if torsion_name not in multiplicity_bitstrings.keys():
                multiplicity_bitstrings[torsion_name] = 0

            for m in multiplicities:
                name = p[0] + '_' + p[1] + '_' + p[2] + '_' + p[3] + '_' + str(m) + '_K'
                if not self.sample_phase:
                    k = pymc.Uniform(name, lower=-20, upper=20, value=0)
                else:
                    k = pymc.Uniform(name, lower=0, upper=20, value=0)

                for i in range(len(param.dihedral_types[p])):
                    if param.dihedral_types[p][i].per == m:
                        multiplicity_bitstrings[torsion_name] += 2 ** (m - 1)
                        if not self.sample_phase:
                            k = pymc.Uniform(name, lower=-20, upper=20, value=param.dihedral_types[p][i].phi_k)
                        else:
                            k = pymc.Uniform(name, lower=0, upper=20, value=param.dihedral_types[p][i].phi_k)
                        break

                self.pymc_parameters[name] = k

                if self.sample_phase:
                    name = p[0] + '_' + p[1] + '_' + p[2] + '_' + p[3] + '_' + str(m) + '_Phase'
                    for i in range(len(param.dihedral_types[p])):
                        if param.dihedral_types[p][i].per == m:
                            if self.continuous_phase:
                                phase = pymc.Uniform(name, lower=0, upper=180.0, value=param.dihedral_types[p][i].phase)
                            else:
                                if param.dihedral_types[p][i].phase == 0:
                                    phase = pymc.DiscreteUniform(name, lower=0, upper=1, value=0)
                                    break

                                if param.dihedral_types[p][i].phase == 180.0:
                                    phase = pymc.DiscreteUniform(name, lower=0, upper=1, value=1)
                                    break
                        else:
                            if self.continuous_phase:
                                phase = pymc.Uniform(name, lower=0, upper=180.0, value=0)
                            else:
                                phase = pymc.DiscreteUniform(name, lower=0, upper=1, value=0)

                    self.pymc_parameters[name] = phase

        if self.rj:
            for torsion_name in multiplicity_bitstrings.keys():
                name = torsion_name + '_multiplicity_bitstring'
                bitstring = pymc.DiscreteUniform(name, lower=0, upper=63, value=multiplicity_bitstrings[torsion_name])
                self.pymc_parameters[name] = bitstring

        self.pymc_parameters['log_sigma'] = pymc.Uniform('log_sigma', lower=-10, upper=3, value=np.log(0.01))
        self.pymc_parameters['sigma'] = pymc.Lambda('sigma',
                                                    lambda log_sigma=self.pymc_parameters['log_sigma']: np.exp(
                                                        log_sigma))
        self.pymc_parameters['precision'] = pymc.Lambda('precision',
                                                        lambda log_sigma=self.pymc_parameters['log_sigma']: np.exp(
                                                            -2 * log_sigma))

        # add missing multiplicity terms to parameterSet so that the system has the same number of parameters
        par.add_missing(self.parameters_to_optimize, param, sample_n5=self.sample_n5)

        @pymc.deterministic
        def mm_energy(pymc_parameters=self.pymc_parameters, param=param):
            mm = np.ndarray(0)
            par.update_param_from_sample(self.parameters_to_optimize, param, model=self, rj=self.rj,
                                         phase=self.sample_phase, n_5=self.sample_n5, continuous=self.continuous_phase)
            for mol in self.frags:
                mol.compute_energy(param, offset=self.pymc_parameters['%s_offset' % mol.topology._residues[0]],
                                   platform=self.platform)
                mm = np.append(mm, mol.mm_energy / kilojoules_per_mole)
            return mm

        size = sum([len(i.qm_energy) for i in self.frags])
        qm_energy = np.ndarray(0)
        for i in range(len(frags)):
             qm_energy = np.append(qm_energy, frags[i].qm_energy)
        #diff_energy = np.ndarray(0)
        #for i in range(len(frags)):
        #    diff_energy = np.append(diff_energy, frags[i].delta_energy)
        self.pymc_parameters['mm_energy'] = mm_energy
        self.pymc_parameters['qm_fit'] = pymc.Normal('qm_fit', mu=self.pymc_parameters['mm_energy'],
                                                     tau=self.pymc_parameters['precision'], size=size, observed=True,
                                                     value=qm_energy)
