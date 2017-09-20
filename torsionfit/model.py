__author__ = 'Chaya D. Stern'

import pymc
import numpy as np
from simtk.unit import kilojoules_per_mole
import torsionfit.database.qmdatabase as TorsionScan
import torsionfit.parameters as par
import warnings
from torsionfit.utils import logger
from collections import OrderedDict
import itertools


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
    def __init__(self, param, frags, stream=None,  platform=None, param_to_opt=None, rj=False, init_random=True, tau='mult'):
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
        init_random: bool
            Randomize starting condition. Default is True. If false, will resort to whatever value is in the parameter set.
        tau: string.
            options are mult or single. When mult, every element in K_m will have its own tau, when single, each K_m will
             have one tau.

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
        if param_to_opt:
            self.parameters_to_optimize = param_to_opt
        else:
            self.parameters_to_optimize = TorsionScan.to_optimize(param, stream)

        multiplicity_bitstrings = dict()

        # offset
        for frag in self.frags:
            name = '%s_offset' % frag.topology._residues[0]
            offset = pymc.Uniform(name, lower=-50, upper=50, value=0)
            self.pymc_parameters[name] = offset

        if tau=='mult':
            value = np.log(np.ones(6)*0.01)
        else:
            value = np.log(0.01)
        for p in self.parameters_to_optimize:
            torsion_name = p[0] + '_' + p[1] + '_' + p[2] + '_' + p[3]

            self.pymc_parameters['log_sigma_k_{}'.format(torsion_name)] = pymc.Uniform('log_sigma_k_{}'.format(torsion_name),
                                                                                       lower=-4.6052, upper=3.453,
                                                                                       value=value)
            self.pymc_parameters['sigma_k_{}'.format(torsion_name)] = pymc.Lambda('sigma_k_{}'.format(torsion_name),
                                                     lambda log_sigma_k=self.pymc_parameters['log_sigma_k_{}'.format(torsion_name)]: np.exp(
                                                        log_sigma_k))
            self.pymc_parameters['precision_k_{}'.format(torsion_name)] = pymc.Lambda('precision_k_{}'.format(torsion_name),
                                lambda log_sigma_k=self.pymc_parameters['log_sigma_k_{}'.format(torsion_name)]: np.exp(
                                                            -2 * log_sigma_k))

            self.pymc_parameters['{}_K'.format(torsion_name)] = pymc.Normal('{}_K'.format(torsion_name), value=np.zeros(6), mu=0,
                                                                            tau=self.pymc_parameters['precision_k_{}'.format(torsion_name)])

            if torsion_name not in multiplicity_bitstrings.keys():
                multiplicity_bitstrings[torsion_name] = 0

        if self.rj:
            for torsion_name in multiplicity_bitstrings.keys():
                name = torsion_name + '_multiplicity_bitstring'
                bitstring = pymc.DiscreteUniform(name, lower=0, upper=63, value=multiplicity_bitstrings[torsion_name])
                self.pymc_parameters[name] = bitstring

        if init_random:
            # randomize initial value
            for parameter in self.pymc_parameters:
                if type(self.pymc_parameters[parameter]) != pymc.CommonDeterministics.Lambda and parameter[:11] != 'log_sigma_k':
                    self.pymc_parameters[parameter].random()
                    logger().info('initial value for {} is {}'.format(parameter, self.pymc_parameters[parameter].value))

        self.pymc_parameters['log_sigma'] = pymc.Uniform('log_sigma', lower=-10, upper=3, value=np.log(0.01))
        self.pymc_parameters['sigma'] = pymc.Lambda('sigma',
                                                    lambda log_sigma=self.pymc_parameters['log_sigma']: np.exp(
                                                        log_sigma))
        self.pymc_parameters['precision'] = pymc.Lambda('precision',
                                                        lambda log_sigma=self.pymc_parameters['log_sigma']: np.exp(
                                                            -2 * log_sigma))

        # Precalculate phis
        n = np.array([1., 2., 3., 4., 5., 6.])
        self.bitmasks = []
        for i in itertools.product((0,1), repeat=6):
            self.bitmasks.append(i)

        cosines = []
        for i, frag in enumerate(frags):
            cosines.append(OrderedDict())
            for t in frag.phis:
                cosines[i][t] = (1 + np.cos(frag.phis[t][:, np.newaxis]*n[:, np.newaxis])).sum(-1)
        self.cosines = cosines

        @pymc.deterministic
        def torsion_energy(pymc_parameters=self.pymc_parameters):
            mm = np.ndarray(0)

            for i, mol in enumerate(self.frags):
                Fourier_sum = np.zeros((mol.n_frames))
                for t in cosines[i]:
                    name = t[0] + '_' + t[1] + '_' + t[2] + '_' + t[3]
                    if self.rj:
                        K = pymc_parameters['{}_K'.format(name)] * self.bitmasks[pymc_parameters['{}_multiplicity_bitstring'.format(name)]]
                    else:
                        K = pymc_parameters['{}_K'.format(name)]
                    Fourier_sum += (K*cosines[i][t]).sum(1)
                Fourier_sum_rel = Fourier_sum - min(Fourier_sum)
                Fourier_sum_rel += pymc_parameters['{}_offset'.format(mol.topology._residues[0])]
                mm = np.append(mm, Fourier_sum)
            return mm

        size = sum([len(i.qm_energy) for i in self.frags])
        residual_energy = np.ndarray(0)
        for i in range(len(frags)):
            residual_energy = np.append(residual_energy, frags[i].delta_energy)

        self.pymc_parameters['torsion_energy'] = torsion_energy
        self.pymc_parameters['qm_fit'] = pymc.Normal('qm_fit', mu=self.pymc_parameters['torsion_energy'],
                                                     tau=self.pymc_parameters['precision'], size=size, observed=True,
                                                     value=residual_energy)
