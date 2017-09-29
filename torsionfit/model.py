__author__ = 'Chaya D. Stern'

import pymc
import numpy as np
import torsionfit.database.qmdatabase as TorsionScan
from torsionfit.utils import logger
from collections import OrderedDict
import itertools


class TorsionFitModel(object):
    """pymc model for sampling torsion parameters.
    This model provides the option to use reversible jump to sample over multiplicity terms

    Attributes:
    ----------
    pymc_parameters: dict() of pymc parameters
    fags: list of QMDatabasees for fragments to fit torsions to
    rj: bool. If True, model uses reversible jump to sample over multiplicity terms. If false, all terms are sampled.
    parameters_to_optimize: list of tuples (dihedrals to optimize)
    models: list of models to sample over.
    inner_sum: list of precalculated inner sum. This is also the gradient.

    """
    def __init__(self, param, frags, stream=None,  param_to_opt=None, rj=False, init_random=True, tau='mult'):
        """

        Parameters
        ----------
        param : Parmed CharmmParameterSet
        frags : list of torsionfit.QMDataBase
        stream : str
            Path to CHARMM stream file. Default None. If None, param_to_opt list must be given. When a stream file is
            specified, param_to_opt is generated if the penalty of the parameters are greater than a threshold.
        param_to_opt : list of tuples of torsions.
            Default None.
        rj : bool
            If True, will use reversible jump to sample Fourier terms. If False, will sample all Ks. Default False
        init_random: bool
            Randomize starting condition. Default is True. If false, will resort to whatever value is in the parameter set.
            Default True
        tau: string.
            options are 'mult' or 'single'. When 'mult', every element in K_m will have its own 'tau', when 'single',
            each K_m will have one tau.
            Default 'mult'

        Returns
        -------
        pymc model

        """

        if type(frags) != list:
            frags = [frags]

        self.pymc_parameters = dict()
        self.frags = frags
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
        elif tau == 'single':
            value = np.log(0.01)
        else:
            raise Exception("Only 'mult' and 'single' are allowed options for tau")

        for p in self.parameters_to_optimize:
            torsion_name = p[0] + '_' + p[1] + '_' + p[2] + '_' + p[3]

            # lower and upper for this distribution are based on empirical data that below this amount the prior is too
            # biased and above the moves are usually rejected.
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
        self.models = []
        for i in itertools.product((0, 1), repeat=6):
            self.models.append(i)

        inner_sum = []
        for i, frag in enumerate(frags):
            inner_sum.append(OrderedDict())
            for t in frag.phis:
                inner_sum[i][t] = (1 + np.cos(frag.phis[t][:, np.newaxis]*n[:, np.newaxis])).sum(-1)
        self.inner_sum = inner_sum

        @pymc.deterministic
        def torsion_energy(pymc_parameters=self.pymc_parameters):
            mm = np.ndarray(0)

            for i, mol in enumerate(self.frags):
                Fourier_sum = np.zeros((mol.n_frames))
                for t in inner_sum[i]:
                    name = t[0] + '_' + t[1] + '_' + t[2] + '_' + t[3]
                    if self.rj:
                        K = pymc_parameters['{}_K'.format(name)] * self.models[pymc_parameters['{}_multiplicity_bitstring'.format(name)]]
                    else:
                        K = pymc_parameters['{}_K'.format(name)]
                    Fourier_sum += (K*inner_sum[i][t]).sum(1)
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
