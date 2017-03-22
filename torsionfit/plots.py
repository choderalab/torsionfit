"""
Plotting module and data exploration for torsionfit

This module contains functions to simplify exploration of torsionfit output in addition to general plotting functions
"""

__author__ = 'Chaya D. Stern'

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pymbar

# global parameter
multiplicities = (1, 2, 3, 4, 6)


def get_parameter_names(model, db, n5=False):
    """
    returns a dictionary that maps torsion name all associated parameters for convenient trace access in pymc.database
    :param model: torsionfit.TorsionFitModel
    :param db: pymc.database (can also be pymc.sampler)

    :return: dictionary mapping torsion name to all associated parameters
    """
    if n5:
        multiplicities = tuple(range(1, 7))
    else:
        multiplicities = (1, 2, 3, 4, 6)
    torsion_parameters = {}
    torsions = model.parameters_to_optimize
    for name in torsions:
        torsion_name = name[0] + '_' + name[1] + '_' + name[2] + '_' + name[3]
        torsion_parameters[torsion_name] = []
        multiplicity_bitstring = torsion_name + '_multiplicity_bitstring'
        torsion_parameters[torsion_name].append(multiplicity_bitstring)
        for m in multiplicities:
            k = torsion_name + '_' + str(m) + '_K'
            torsion_parameters[torsion_name].append(k)
            phase = torsion_name + '_' + str(m) + '_Phase'
            torsion_parameters[torsion_name].append(phase)
    return torsion_parameters


def get_multiplicity_traces(torsion_parameters, db, n5=False):
    """
    returns traces for the multiplicity terms for all torsions in (0, 1)

    :param torsion_parameters: dict mapping torsion name to parameters for that torsion or name of torsion
    :param db: pymc.database

    :return: dict mapping torsion name to multiplicity terms trace
    """
    if n5:
        multiplicities = tuple(range(1, 7))
    else:
        multiplicities = (1, 2, 3, 4, 6)
    if type(torsion_parameters) == str:
        torsion_parameters = [torsion_parameters]
    else:
        torsion_parameters = torsion_parameters.keys()
    multiplicity_traces = {}
    for torsion_name in torsion_parameters:
        multiplicity_bitstring = torsion_name + '_multiplicity_bitstring'
        for m in multiplicities:
            multiplicity_traces[torsion_name + '_' + str(m)] = []
            for i in db.trace(multiplicity_bitstring)[:]:
                if 2**(m-1) & int(i):
                    multiplicity_traces[torsion_name + '_' + str(m)].append(1)
                else:
                    multiplicity_traces[torsion_name + '_' + str(m)].append(0)
    return multiplicity_traces


def get_statistics(db, torsion_parameters):
    """
    uses pymbar.timeseries.detectEquilibration module to get equilibration time, statistical inefficiency and effective
    samples for each trace. Returns a dictionary that maps all parameters to statistics.

    :param db: pymc.database (can also use pymc.sampler)
    :param torsion_parameters: dict mapping torsion name to associated parameters

    :return: dict that maps parameters to statistics
    """

    statistics = {}
    for parameters in torsion_parameters:
        for param in torsion_parameters[parameters]:
            statistics[param] = pymbar.timeseries.detectEquilibration(db.trace(param)[:])

    return statistics


def trace_plots(name, db, markersize, statistics=False, multiplicity_traces=False, continuous=False, filename=None):
    """
    Generate trace plot for all parameters of a given torsion

    :param name: str. name of torsion parameter A_B_C_D where A, B, C, and D are atom types.
    :param db: pymc.database (can also use pymc.sampler)
    :param markersize: int.
    :param statistics: dict that maps parameters to statistics from pymbar.timeseries.detectEquilibrium. Default: False
    :param multiplicity_traces: dict that maps multiplicity term to (0,1) trace. Default is False.
    """

    if not multiplicity_traces:
        try:
            multiplicity_traces = get_multiplicity_traces(torsion_parameters=name, db=db)
        except KeyError:
            pass

    pp = PdfPages('%s_traces.pdf' % name)
    fig = plt.figure()

    axes_k = plt.subplot(9, 2, 1)
    plt.plot(db.trace(name + '_' + str(1) + '_K')[:], 'k.', markersize=markersize, label='K')
    plt.title(name, fontweight='bold')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '1' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(1) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(0, 20)
    plt.ylabel('kJ/mole')
    plt.xticks([])
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.yticks([0, 20])

    axes_phase = plt.subplot(9, 2, 3)
    plt.plot(db.trace(name + '_' + str(1) + '_Phase')[:], '.', markersize=markersize, label='Phase')
    if continuous:
        plt.ylim(-1.0, 181)
        plt.yticks([1, 180])
    else:
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])

    axes_n = plt.subplot(9, 2, 5)
    try:
        plt.plot(multiplicity_traces[name + '_' + str(1)], 'k.', markersize=markersize, label='1')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])
        plt.xticks([])
    except:
        pass

    axes_k = plt.subplot(9, 2, 7)
    plt.plot(db.trace(name + '_' + str(2) + '_K')[:], 'k.', markersize=markersize, label='K')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '2' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(2) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(0, 20)
    plt.ylabel('kJ/mole')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])
    plt.yticks([0, 20])

    axes_phase = plt.subplot(9, 2, 9)
    plt.plot(db.trace(name + '_' + str(2) + '_Phase')[:], '.', markersize=markersize, label='Phase')
    if continuous:
        plt.ylim(-1.0, 181)
        plt.yticks([1, 180])
    else:
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])

    axes_n = plt.subplot(9, 2, 11)
    try:
        plt.plot(multiplicity_traces[name + '_' + str(2)], 'k.', markersize=markersize, label='2')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])
        plt.xticks([])
    except:
        pass

    axes_k = plt.subplot(9, 2, 13)
    plt.plot(db.trace(name + '_' + str(3) + '_K')[:], 'k.', markersize=markersize, label='K')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '3' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(3) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(0, 20)
    plt.ylabel('kJ/mole')
    plt.xticks([])
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.yticks([0, 20])

    axes_phase = plt.subplot(9, 2, 15)
    plt.plot(db.trace(name + '_' + str(3) + '_Phase')[:], '.', markersize=markersize, label='Phase')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])
    if continuous:
        plt.ylim(-1.0, 181)
        plt.yticks([1, 180])
    else:
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])

    axes_n = plt.subplot(9, 2, 17)
    try:
        plt.plot(multiplicity_traces[name + '_' + str(3)], 'k.', markersize=markersize, label='3')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])
        plt.xlabel('mcmc steps')
    except:
        pass

    axes_k = plt.subplot(9, 2, 2)
    plt.title(name, fontweight='bold')
    plt.plot(db.trace(name + '_' + str(4) + '_K')[:], 'k.', markersize=markersize, label='K')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '4' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(4) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(0, 20)
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])
    plt.yticks([])

    axes_phase = plt.subplot(9, 2, 4)
    plt.plot(db.trace(name + '_' + str(4) + '_Phase')[:], '.', markersize=markersize, label='Phase')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    if continuous:
        plt.ylim(-1.0, 181)
        plt.yticks([1, 180])
    else:
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])
    plt.xticks([])
    plt.yticks([])

    try:
        axes_n = plt.subplot(9, 2, 6)
        plt.plot(multiplicity_traces[name + '_' + str(4)], 'k.', markersize=markersize, label='4')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([])
        plt.xticks([])
    except:
        pass

    axes_k = plt.subplot(9, 2, 8)
    plt.plot(db.trace(name + '_' + str(6) + '_K')[:], 'k.', markersize=markersize, label='K')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '6' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(6) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(0, 20)
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])
    plt.yticks([])

    axes_phase = plt.subplot(9, 2, 10)
    plt.plot(db.trace(name + '_' + str(6) + '_Phase')[:], '.', markersize=markersize, label='Phase')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    if continuous:
        plt.ylim(-1.0, 181)
        plt.yticks([1, 180])
    else:
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])
    plt.xticks([])
    plt.yticks([])

    axes_n = plt.subplot(9, 2, 12)
    try:
        plt.plot(multiplicity_traces[name + '_' + str(6)], 'k.', markersize=markersize, label='6')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([])
        plt.xlabel('mcmc steps')
    except:
        pass

    if filename is None:
        fig.savefig('%s_traces.pdf' % name)
    else:
        fig.savefig(filename)
    pp.savefig(fig, dpi=80)
    pp.close()


def trace_no_phase(name, db, markersize, statistics=False, multiplicity_traces=False, ymin=-20, ymax=20, filename=None):
    """
    Generate trace plot for all parameters of a given torsion

    :param name: str. name of torsion parameter A_B_C_D where A, B, C, and D are atom types.
    :param db: pymc.database (can also use pymc.sampler)
    :param markersize: int.
    :param statistics: dict that maps parameters to statistics from pymbar.timeseries.detectEquilibrium. Default: False
    :param multiplicity_traces: dict that maps multiplicity term to (0,1) trace. Default is False.
    """

    if not multiplicity_traces:
        try:
            multiplicity_traces = get_multiplicity_traces(torsion_parameters=name, db=db)
        except KeyError:
            pass

    pp = PdfPages('%s_traces.pdf' % name)
    fig = plt.figure()

    axes_k = plt.subplot(5, 2, 1)
    plt.plot(db.trace(name + '_' + str(1) + '_K')[:], 'k.', markersize=markersize, label='K')
    plt.title(name, fontweight='bold')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '1' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(1) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.ylabel('kJ/mole')
    plt.xticks([])
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.yticks([ymin, 0, ymax])

    axes_n = plt.subplot(5, 2, 2)
    try:
        plt.plot(multiplicity_traces[name + '_' + str(1)], 'k.', markersize=markersize, label='1')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])
        plt.xticks([])
    except:
        pass

    axes_k = plt.subplot(5, 2, 3)
    plt.plot(db.trace(name + '_' + str(2) + '_K')[:], 'k.', markersize=markersize, label='K')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '2' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(2) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.ylabel('kJ/mole')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])
    plt.yticks([ymin, 0, ymax])

    axes_n = plt.subplot(5, 2, 4)
    try:
        plt.plot(multiplicity_traces[name + '_' + str(2)], 'k.', markersize=markersize, label='2')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([0, 1])
        plt.xticks([])
    except:
        pass

    axes_k = plt.subplot(5, 2, 5)
    plt.plot(db.trace(name + '_' + str(3) + '_K')[:], 'k.', markersize=markersize, label='K')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '3' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(3) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.ylabel('kJ/mole')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.yticks([ymin, 0,  ymax])
    plt.xlabel('mcmc steps')


    axes_n = plt.subplot(5, 2, 6)
    plt.title(name, fontweight='bold')
    try:
        plt.plot(multiplicity_traces[name + '_' + str(3)], 'k.', markersize=markersize, label='3')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([])
        plt.xticks([])
    except:
        pass

    axes_k = plt.subplot(5, 2, 7)
    plt.plot(db.trace(name + '_' + str(4) + '_K')[:], 'k.', markersize=markersize, label='K')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '4' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(4) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])
    plt.yticks([ymin, 0, ymax])

    axes_n = plt.subplot(5, 2, 8)
    try:
        plt.plot(multiplicity_traces[name + '_' + str(4)], 'k.', markersize=markersize, label='4')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([])
        plt.xticks([])
    except:
        pass

    axes_k = plt.subplot(5, 2, 9)
    plt.plot(db.trace(name + '_' + str(6) + '_K')[:], 'k.', markersize=markersize, label='K')
    if statistics:
        axes_k.axvline(statistics[name + '_' + '6' + '_K'][0], color='red', lw=1)
    else:
        axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(6) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.yticks([])
    plt.xticks([ymin, 0, ymax])

    axes_n = plt.subplot(5, 2, 10)
    try:
        plt.plot(multiplicity_traces[name + '_' + str(6)], 'k.', markersize=markersize, label='6')
        plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
        plt.ylim(-0.1, 1.1)
        plt.yticks([])
        plt.xlabel('mcmc steps')
    except:
        pass

    if not filename:
        fig.savefig('%s_traces.pdf' % name)
    else:
        fig.savefig(filename)
    pp.savefig(fig, dpi=80)
    pp.close()

def trace_no_phase_n5(name, db, markersize, statistics=False, equil=True,  multiplicity_traces=False, ymin=-20, ymax=20, filename=None):
    """
    Generate trace plot for all parameters of a given torsion

    :param name: str. name of torsion parameter A_B_C_D where A, B, C, and D are atom types.
    :param db: pymc.database (can also use pymc.sampler)
    :param markersize: int.
    :param statistics: dict that maps parameters to statistics from pymbar.timeseries.detectEquilibrium. Default: False
    :param multiplicity_traces: dict that maps multiplicity term to (0,1) trace. Default is False.
    """

    if not multiplicity_traces:
        multiplicity_traces = get_multiplicity_traces(torsion_parameters=name, db=db, n5=True)

    pp = PdfPages('%s_traces.pdf' % name)
    fig = plt.figure()

    axes_k = plt.subplot(6, 2, 1)
    plt.plot(db.trace(name + '_' + str(1) + '_K')[:], 'k.', markersize=markersize, label='K')
    plt.title(name, fontweight='bold')
    if equil:
        if statistics:
            axes_k.axvline(statistics[name + '_' + '1' + '_K'][0], color='red', lw=1)
        else:
            axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(1) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.ylabel('kJ/mole')
    plt.xticks([])
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.yticks([ymin, 0, ymax])

    axes_n = plt.subplot(6, 2, 2)
    plt.title(name, fontweight='bold')
    plt.plot(multiplicity_traces[name + '_' + str(1)], 'k.', markersize=markersize, label='1')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.ylim(-0.1, 1.1)
    plt.yticks([0, 1])
    plt.xticks([])

    axes_k = plt.subplot(6, 2, 3)
    plt.plot(db.trace(name + '_' + str(2) + '_K')[:], 'k.', markersize=markersize, label='K')
    if equil:
        if statistics:
            axes_k.axvline(statistics[name + '_' + '2' + '_K'][0], color='red', lw=1)
        else:
            axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(2) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.ylabel('kJ/mole')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])
    plt.yticks([ymin, 0, ymax])

    axes_n = plt.subplot(6, 2, 4)
    plt.plot(multiplicity_traces[name + '_' + str(2)], 'k.', markersize=markersize, label='2')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.ylim(-0.1, 1.1)
    plt.yticks([0, 1])
    plt.xticks([])

    axes_k = plt.subplot(6, 2, 5)
    plt.plot(db.trace(name + '_' + str(3) + '_K')[:], 'k.', markersize=markersize, label='K')
    if equil:
        if statistics:
            axes_k.axvline(statistics[name + '_' + '3' + '_K'][0], color='red', lw=1)
        else:
            axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(3) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.ylabel('kJ/mole')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.yticks([ymin, 0,  ymax])
    plt.xticks([])

    axes_n = plt.subplot(6, 2, 6)
    plt.plot(multiplicity_traces[name + '_' + str(3)], 'k.', markersize=markersize, label='3')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.ylim(-0.1, 1.1)
    plt.yticks([0, 1])
    plt.xticks([])

    axes_k = plt.subplot(6, 2, 7)
    plt.plot(db.trace(name + '_' + str(4) + '_K')[:], 'k.', markersize=markersize, label='K')
    if equil:
        if statistics:
            axes_k.axvline(statistics[name + '_' + '4' + '_K'][0], color='red', lw=1)
        else:
            axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(4) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])
    plt.yticks([ymin, 0, ymax])
    plt.ylabel('KJ/mol')

    axes_n = plt.subplot(6, 2, 8)
    plt.plot(multiplicity_traces[name + '_' + str(4)], 'k.', markersize=markersize, label='4')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.ylim(-0.1, 1.1)
    plt.yticks([0, 1])
    plt.xticks([])

    axes_k = plt.subplot(6, 2, 9)
    plt.plot(db.trace(name + '_' + str(5) + '_K')[:], 'k.', markersize=markersize, label='K')
    if equil:
        if statistics:
            axes_k.axvline(statistics[name + '_' + '5' + '_K'][0], color='red', lw=1)
        else:
            axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(5) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xticks([])
    plt.yticks([ymin, 0, ymax])
    plt.ylabel('KJ/mol')

    axes_n = plt.subplot(6, 2, 10)
    plt.plot(multiplicity_traces[name + '_' + str(5)], 'k.', markersize=markersize, label='5')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.ylim(-0.1, 1.1)
    plt.yticks([0, 1])
    plt.xticks([])

    axes_k = plt.subplot(6, 2, 11)
    plt.plot(db.trace(name + '_' + str(6) + '_K')[:], 'k.', markersize=markersize, label='K')
    if equil:
        if statistics:
            axes_k.axvline(statistics[name + '_' + '6' + '_K'][0], color='red', lw=1)
        else:
            axes_k.axvline(pymbar.timeseries.detectEquilibration(db.trace(name + '_' + str(6) + '_K')[:])[0], color='red',
                       lw=1)
    plt.ylim(ymin, ymax)
    plt.yticks([ymin, 0, ymax])
    plt.ylabel('KJ/mol')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.xlabel('mcmc steps')

    axes_n = plt.subplot(6, 2, 12)
    plt.plot(multiplicity_traces[name + '_' + str(6)], 'k.', markersize=markersize, label='6')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.ylim(-0.1, 1.1)
    plt.yticks([0, 1])
    plt.xlabel('mcmc steps')

    if not filename:
        fig.savefig('%s_traces.pdf' % name)
    else:
        fig.savefig(filename)
    pp.savefig(fig, dpi=80)
    pp.close()

def marg_mult(model, db, samples, burn=0, filename=None, n5=False):
    """
    generates histogram for marginal distribution of posterior multiplicities.

    :param model: TorsionFitModel
    :param db: pymc.database for model
    :param samples: length of trace
    :param burn: int. number of steps to skip
    :param filename: filename for plot to save
    """
    if n5:
        multiplicities = tuple(range(1, 7))
    else:
        multiplicities = (1, 2, 3, 4, 6)
    mult_bitstring = []
    for i in model.pymc_parameters.keys():
        if i.split('_')[-1] == 'bitstring':
            mult_bitstring.append(i)
    if n5:
        histogram = np.zeros((len(mult_bitstring), samples, 5))
    else:
        histogram = np.zeros((len(mult_bitstring), samples, 5))

    for m, torsion in enumerate(mult_bitstring):
        for i, j in enumerate(db.trace('%s' % torsion)[burn:]):
            for k, l in enumerate(multiplicities):
                if 2**(l-1) & int(j):
                    histogram[m][i][k] = 1

    plt.matshow(histogram.sum(1), cmap='jet',  extent=[0, 5, 0, 20]), plt.colorbar()
    plt.yticks([])
    plt.xlabel('multiplicity term')
    plt.ylabel('torsion')
    if filename:
        plt.savefig(filename)
