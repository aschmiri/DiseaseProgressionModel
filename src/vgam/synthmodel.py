'''
A class to provide modelling functionality for synthetic examples.

@author:     Alexander Schmidt-Richberg
@copyright:  2014 Imperial College London. All rights reserved.
@contact:    a.schmidt-richberg@imperial.ac.uk
'''
import math
import random
import numpy as np
from common import log as log
from scipy import stats


class SynthModel(object):
    '''
    A class to simulate synthetic models.
    '''
    PROGRESS_RANGE = 2000

    _models = {}
    # CDR-SB model following [Hennemann 09]
    #   Group CDR = 0.5: 3.2 (0.9)
    #   Group CDR = 1-2: 7.9 (2.6)
    _models.update({'synth_cdrsb': {'shape': 'exp',
                                    'exp_asymptote': 0,
                                    'exp_a': 3.2,
                                    'exp_b': np.log(7.9 / 3.2) / float(365 * 6),
                                    'offset': 365 * -3,
                                    'noise': 'gamma',
                                    'noise_sigma_mci': 0.9,
                                    'noise_sigma_ad': 2.6}})

    # MMSE model following [Hennemann 09]
    #   Group CDR = 0.5: 23.2 (2.6)
    #   Group CDR = 1-2: 18.4 (4.0)
    _models.update({'synth_mmse': {'shape': 'exp',
                                   'exp_asymptote': 30,
                                   'exp_a': 6.8 * -1,
                                   'exp_b': np.log(11.6 / 6.8) / float(365 * 6),
                                   'offset': 365 * -3,
                                   'noise': 'gamma',
                                   'noise_sigma_mci': 2.6,
                                   'noise_sigma_ad': 4.0}})

    # Hippocampus volume model following [Coley 09]
    #   Mean volume: 3653 mm^2, Atrophy rate: 3.3%
    _models.update({'synth_hipp': {'shape': 'pow',
                                   'pow_a': 3693,
                                   'pow_b': 0.965,
                                   'offset': 0.0,
                                   'noise': 'gaussian',
                                   'noise_sigma_mci': 572.0,
                                   'noise_sigma_ad': 572.0}})

    # Whole brain volume model following [Coley 09]
    #   Mean volume: 1487 ml, Atrophy rate: 0.3%
    _models.update({'synth_brain': {'shape': 'pow',
                                    'pow_a': 1480,
                                    'pow_b': 0.986,
                                    'offset': 0.0,
                                    'noise': 'gaussian',
                                    'noise_sigma_mci': 92.0,
                                    'noise_sigma_ad': 92.0}})

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
    def get_random_progress(sampling='uniform', sample_range=None):
        sample_range = SynthModel.PROGRESS_RANGE if sample_range is None else sample_range
        if sampling == 'uniform':
            return int(random.uniform(-sample_range, sample_range))
        elif sampling == 'triangular':
            return int(random.triangular(-sample_range, sample_range, 0))
        elif sampling == 'longitudinal':
            return int(random.uniform(-sample_range, 0))
        else:
            print log.ERROR, 'Unknown sampling: {0}'.format(sampling)
            return None

    ############################################################################
    #
    # get_distributed_value()
    #
    ############################################################################
    @staticmethod
    def get_distributed_value(biomarker, progress, fixed_sigma=None, cdf=None):
        if SynthModel._models[biomarker]['noise'] == 'gaussian':
            median = SynthModel.get_median(biomarker, progress)
            sigma = fixed_sigma if fixed_sigma is not None else SynthModel.get_sigma(biomarker, progress)
            if cdf is None:
                return median + random.gauss(0.0, sigma)
            else:
                return stats.norm.ppf(cdf, loc=median, scale=sigma)
        elif SynthModel._models[biomarker]['noise'] == 'gamma':
            k, theta = SynthModel._get_gamma_parameters(biomarker, progress)
            if cdf is None:
                return SynthModel._transorm_coordinates(biomarker, random.gammavariate(k, theta))
            else:
                return SynthModel._transorm_coordinates(biomarker, stats.gamma.ppf(cdf, k, scale=theta))

    ############################################################################
    #
    # get_median()
    #
    ############################################################################
    @staticmethod
    def get_median(biomarker, progress):
        if SynthModel._models[biomarker]['shape'] == 'lin':
            slope = SynthModel._models[biomarker]['lin_a']
            offset = SynthModel._models[biomarker]['offset']
            return slope * (progress - offset)
        elif SynthModel._models[biomarker]['shape'] == 'exp':
            asymptote = SynthModel._models[biomarker]['exp_asymptote']
            a = SynthModel._models[biomarker]['exp_a']
            b = SynthModel._models[biomarker]['exp_b']
            offset = SynthModel._models[biomarker]['offset']
            return asymptote + a * np.exp(b * (progress - offset))
        elif SynthModel._models[biomarker]['shape'] == 'pow':
            a = SynthModel._models[biomarker]['pow_a']
            b = SynthModel._models[biomarker]['pow_b']
            offset = SynthModel._models[biomarker]['offset']
            return a * b ** ((progress - offset) / 365.0)
        elif SynthModel._models[biomarker]['shape'] == 'sigmoid':
            a = SynthModel._models[biomarker]['sig_a']
            min_val = SynthModel._models[biomarker]['sig_min']
            max_val = SynthModel._models[biomarker]['sig_max']
            offset = SynthModel._models[biomarker]['offset']
            return min_val + (max_val - min_val) / (1.0 + np.exp(a * (progress - offset)))

    ############################################################################
    #
    # get_sigma()
    #
    ############################################################################
    @staticmethod
    def get_sigma(biomarker, progress):
        '''
        Get an interpolated value for sigma as an exponential function, such that
          sigma(-3 * 365) = noise_sigma_mci
          sigma(3 * 365) = noise_sigma_ad
        '''
        sigma_mci = SynthModel._models[biomarker]['noise_sigma_mci']
        sigma_ad = SynthModel._models[biomarker]['noise_sigma_ad']
        offset = SynthModel._models[biomarker]['offset']
        if sigma_mci == sigma_ad:
            return sigma_mci
        else:
            b = np.log(sigma_ad / sigma_mci) / float(6 * 365)
            return sigma_mci * np.exp(b * (progress - offset))

    ############################################################################
    #
    # get_probability()
    #
    ############################################################################
    @staticmethod
    def get_probability(biomarker, progress, value):
        if SynthModel._models[biomarker]['noise'] == 'gaussian':
            sigma = SynthModel.get_sigma(biomarker, progress)
            return stats.norm.pdf(value, scale=sigma, loc=SynthModel.get_median(biomarker, progress))
        elif SynthModel._models[biomarker]['noise'] == 'gamma':
            k, theta = SynthModel._get_gamma_parameters(biomarker, progress)
            value = SynthModel._transorm_coordinates(biomarker, value)
            return stats.gamma.pdf(value, k, scale=theta)

    ############################################################################
    #
    # get_cumulated_probability()
    #
    ############################################################################
    @staticmethod
    def get_cumulated_probability(biomarker, progress, value):
        if SynthModel._models[biomarker]['noise'] == 'gaussian':
            sigma = SynthModel.get_sigma(biomarker, progress)
            return stats.norm.cdf(value, scale=sigma, loc=SynthModel.get_median(biomarker, progress))
        elif SynthModel._models[biomarker]['noise'] == 'gamma':
            k, theta = SynthModel._get_gamma_parameters(biomarker, progress)
            value = SynthModel._transorm_coordinates(biomarker, value)
            return stats.gamma.cdf(value, k, scale=theta)

    ############################################################################
    #
    # get_probability()
    #
    ############################################################################
    @staticmethod
    def _get_gamma_parameters(biomarker, progress):
        '''
        Get the values of k and theta for the given biomarker and progress.
        It is calculated based on median and standard deviation sigma using:
          theta = sqrt(sigma^2/k)
        and
          median = k * theta * (3 * k - 0.8) / (3 * k + 0.2)
        '''
        if SynthModel._models[biomarker]['noise'] is not 'gamma':
            print log.ERROR, 'Noise model is not gamma destributed!'
            return None, None
        else:
            median = SynthModel.get_median(biomarker, progress)
            median = SynthModel._transorm_coordinates(biomarker, median)
            sigma = SynthModel.get_sigma(biomarker, progress)

            # Get k from polynomial roots
            m2 = median * median
            s2 = sigma * sigma
            p = [9 * s2, -9 * m2, -(4.8 * s2 + 1.2 * m2), -0.04 * m2, 0.64 * s2]
            roots = np.roots(p)
            for root in roots:
                r = np.real(root)
                if r > 0 and np.imag(r) == 0.0:
                    k = r
                    break

            # Get theta
            theta = (median / k) * (3 * k + 0.2) / (3 * k - 0.8)

            # Return results
            return k, theta

    ############################################################################
    #
    # get_probability()
    #
    ############################################################################
    @staticmethod
    def _transorm_coordinates(biomarker, value):
        if SynthModel._models[biomarker]['shape'] == 'exp':
            asymptote = SynthModel._models[biomarker]['exp_asymptote']
            a = SynthModel._models[biomarker]['exp_a']
            return asymptote + math.copysign(value, a)
        else:
            return value
