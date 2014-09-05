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
                               'noise_std': 0.5}})
    _models.update({'synth2': {'shape': 'exp',
                               'slope': 0.0005,
                               'offset': 0,
                               'noise': 'gaussian',
                               'noise_std': 1.0}})
    _models.update({'synth3': {'shape': 'exp',
                               'slope': 0.001,
                               'offset': 1000,
                               'noise': 'gaussian',
                               'noise_std': 0.5}})
    _models.update({'synth4': {'shape': 'exp',
                               'slope': 0.001,
                               'offset': 1000,
                               'noise': 'gaussian',
                               'noise_std': 1.0}})

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
    def get_random_progress(uniform_progress=False):
        if uniform_progress:
            return random.randint(-SynthModel.PROGRESS_RANGE, SynthModel.PROGRESS_RANGE)
        else:
            vol = random.random()
            if vol < 0.5:
                return int(-SynthModel.PROGRESS_RANGE * (1 - math.sqrt(2 * vol)))
            else:
                return int(SynthModel.PROGRESS_RANGE * (1 - math.sqrt(2 * (1 - vol))))

    ############################################################################
    #
    # get_noisy_value()
    #
    ############################################################################
    @staticmethod
    def get_noisy_value(biomarker, progress):
        return SynthModel.get_value(biomarker, progress) + SynthModel.get_noise(biomarker)

    ############################################################################
    #
    # get_value()
    #
    ############################################################################
    @staticmethod
    def get_value(biomarker, progress):
        if SynthModel._models[biomarker]['shape'] == 'exp':
            slope = SynthModel._models[biomarker]['slope']
            offset = SynthModel._models[biomarker]['offset']
            return np.exp(slope * (progress - offset))

    ############################################################################
    #
    # get_random_value()
    #
    ############################################################################
    @staticmethod
    def get_noise(biomarker):
        if SynthModel._models[biomarker]['noise'] == 'gaussian':
            std = SynthModel._models[biomarker]['noise_std']
            return random.gauss(0, std)

    ############################################################################
    #
    # get_random_value()
    #
    ############################################################################
    @staticmethod
    def get_probability(biomarker, progress, value):
        if SynthModel._models[biomarker]['noise'] == 'gaussian':
            std = SynthModel._models[biomarker]['noise_std']
            offset = value - SynthModel.get_value(biomarker, progress)
            return np.exp(-0.5 * (offset / std) ** 2) / (std * math.sqrt(2 * math.pi))
