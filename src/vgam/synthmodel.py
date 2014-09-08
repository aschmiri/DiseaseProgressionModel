'''
A class to provide modelling functionality for synthetic examples.

@author:     Alexander Schmidt-Richberg
@copyright:  2014 Imperial College London. All rights reserved.
@contact:    a.schmidt-richberg@imperial.ac.uk
'''
import math
import random
import numpy as np


class SynthModel(object):
    '''
    TODO: classdocs
    '''
    PROGRESS_RANGE = 2000
    _models = {}
    _models.update({'synth1': {'shape': 'exp',
                               'slope': 0.0005,
                               'offset': 0,
                               'noise': 'gaussian',
                               'gaussian_std': 0.5}})
    _models.update({'synth2': {'shape': 'exp',
                               'slope': 0.0005,
                               'offset': 0,
                               'noise': 'gaussian',
                               'gaussian_std': 1.0}})
    _models.update({'synth3': {'shape': 'exp',
                               'slope': 0.0005,
                               'offset': (-1000),
                               'noise': 'gamma',
                               'gamma_alpha': 10}})
    _models.update({'synth4': {'shape': 'lin',
                               'slope': 0.0012,
                               'offset': (-3000),
                               'noise': 'gamma',
                               'gamma_alpha': 10}})

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self):
        pass

    ############################################################################
    #
    # get_biomarker_names()
    #
    ############################################################################
    @staticmethod
    def get_biomarker_names():
        return SynthModel._models.keys()

    ############################################################################
    #
    # get_random_progress()
    #
    ############################################################################
    @staticmethod
    def get_random_progress(uniform_progression=False):
        if uniform_progression:
            return int(random.uniform(-SynthModel.PROGRESS_RANGE, SynthModel.PROGRESS_RANGE))
        else:
            return int(random.triangular(-SynthModel.PROGRESS_RANGE, SynthModel.PROGRESS_RANGE, 0))

    ############################################################################
    #
    # get_distributed_value()
    #
    ############################################################################
    @staticmethod
    def get_distributed_value(biomarker, progress):
        if SynthModel._models[biomarker]['noise'] == 'gaussian':
            median = SynthModel.get_median(biomarker, progress)
            std = SynthModel._models[biomarker]['gaussian_std']
            return median + random.gauss(0, std)
        elif SynthModel._models[biomarker]['noise'] == 'gamma':
            alpha = SynthModel._models[biomarker]['gamma_alpha']
            median = SynthModel.get_median(biomarker, progress)
            beta = median / alpha * (3 * alpha + 0.2) / (3 * alpha - 0.8)
            return random.gammavariate(alpha, beta)

    ############################################################################
    #
    # get_median()
    #
    ############################################################################
    @staticmethod
    def get_median(biomarker, progress):
        if SynthModel._models[biomarker]['shape'] == 'lin':
            slope = SynthModel._models[biomarker]['slope']
            offset = SynthModel._models[biomarker]['offset']
            return slope * (progress - offset)
        elif SynthModel._models[biomarker]['shape'] == 'exp':
            slope = SynthModel._models[biomarker]['slope']
            offset = SynthModel._models[biomarker]['offset']
            return np.exp(slope * (progress - offset))

    ############################################################################
    #
    # get_probability()
    #
    ############################################################################
    @staticmethod
    def get_probability(biomarker, progress, value):
        if SynthModel._models[biomarker]['noise'] == 'gaussian':
            x = value - SynthModel.get_median(biomarker, progress)
            std = SynthModel._models[biomarker]['gaussian_std']
            return np.exp(-0.5 * (x / std) ** 2) / (std * math.sqrt(2 * math.pi))
        elif SynthModel._models[biomarker]['noise'] == 'gamma':
            if value < 0:
                return 0
            else:
                alpha = SynthModel._models[biomarker]['gamma_alpha']
                median = SynthModel.get_median(biomarker, progress)
                beta = median / alpha * (3 * alpha + 0.2) / (3 * alpha - 0.8)
                return (value ** (alpha - 1) * math.exp(-value / beta)) / (math.gamma(alpha) * beta ** alpha)
