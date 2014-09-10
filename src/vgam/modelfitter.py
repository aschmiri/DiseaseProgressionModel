'''
A class to provide progression modelling functionality.

@author:     Alexander Schmidt-Richberg
@copyright:  2014 Imperial College London. All rights reserved.
@contact:    a.schmidt-richberg@imperial.ac.uk
'''
import numpy as np
from common import log as log
from common import adni_tools as adni
from vgam.progressionmodel import ProgressionModel
from vgam.progressionmodel import MultiBiomarkerProgressionModel

TEST_SCALE_MIN = 0.1
TEST_SCALE_MAX = 3.0
TEST_SCALE_STEP = 0.1

TEST_DPI_MIN = -2500
TEST_DPI_MAX = 2500
TEST_DPI_STEP = 10

TEST_DPR_MIN = 0.1
TEST_DPR_MAX = 3.0
TEST_DPR_STEP = 0.1


class ModelFitter(object):
    '''
    TODO: classdocs
    '''

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self, model):
        '''
        Initialise the progression model fitter.
        '''
        if isinstance(model, ProgressionModel):
            self.model = model
            print log.INFO, 'Initialised model fitter with {0}'.format(model.biomarker)
        elif isinstance(model, MultiBiomarkerProgressionModel):
            self.model = model
            print log.INFO, 'Initialised model fitter with {0}'.format(model.models.keys())
        else:
            print log.ERROR, 'Invalid model type!'

    ############################################################################
    #
    # get_scaled_measurements()
    #
    ############################################################################
    def get_scaled_measurements(self, measurements):
        ''' Return the optimally scaled measurements.

        Arguments:
        measurements --

        Returns:
        The optimally scaled measurements
        '''
        for rid in measurements:
            print log.INFO, 'Estimating optimal scaling for subject {0}...'.format(rid)

            # Collect samples
            samples = {}
            for viscode in measurements[rid]:
                samples.update({viscode: measurements[rid][viscode]})

            # Get and save scaling
            scaling = self.__get_scaling_for_samples(samples)
            print log.RESULT, 'Optimal scaling value:', scaling

            # Update all progresses
            for viscode in measurements[rid]:
                progress = adni.safe_cast(measurements[rid][viscode]['progress'], int)
                measurements[rid][viscode].update({'scaling': scaling})
                measurements[rid][viscode].update({'progress': progress * scaling})

        return measurements

    ############################################################################
    #
    # __get_scaling_for_samples()
    #
    ############################################################################
    def __get_scaling_for_samples(self, samples):
        ''' Return the optimal scaling value for a subject given a number of samples
        and a set of biomarkers.

        Arguments:
        samples --

        Returns:
        The estimated DPI and DPR
        '''
        test_scalings = np.arange(TEST_SCALE_MIN, TEST_SCALE_MAX, TEST_SCALE_STEP)

        # Compute probability for each scaling
        probs = []
        for scaling in test_scalings:
            prob = 1.0
            for viscode in samples:
                scaled_progress = samples[viscode]['progress'] * scaling
                prob *= self.model.get_probability_value(samples[viscode], scaled_progress)
            # prob /= len(samples)
            probs.append(prob)

        # Sanity check
        if np.sum(probs) == 0.0:
            print log.WARNING, 'All probabilities equal zero!'
            return None

        # Find the scaling with the highest probability
        return test_scalings[np.argmax(probs)]

    ############################################################################
    #
    # get_dpi_for_samples()
    #
    ############################################################################
    def get_dpi_for_samples(self, samples):
        ''' Return the estimated DPI of a subject given a number of samples and
        a set of biomarkers.

        Arguments:
        samples --

        Returns:
        The estimated DPI
        '''
        test_dpis = np.arange(TEST_DPI_MIN, TEST_DPI_MAX, TEST_DPI_STEP)

        # Compute probability for each scaling
        probs = []
        for dpi in test_dpis:
            prob = 1.0
            for viscode in samples:
                offset = samples[viscode]['scantime']
                prob *= self.model.get_probability_value(samples[viscode], dpi + offset)
            # prob /= len(samples)
            probs.append(prob)

        # Sanity check
        if np.sum(probs) == 0.0:
            print log.WARNING, 'All probabilities equal zero!'
            return None

        # Find the DPI with the highest probability
        return test_dpis[np.argmax(probs)]

    ############################################################################
    #
    # get_dpi_dpr_for_samples()
    #
    ############################################################################
    def get_dpi_dpr_for_samples(self, samples):
        ''' Return the estimated DPI and DPR of a subject given a number of samples
        and a set of biomarkers.

        Arguments:
        samples --

        Returns:
        The estimated DPI and DPR
        '''
        test_dprs = np.arange(TEST_DPR_MIN, TEST_DPR_MAX, TEST_DPR_STEP)
        test_dpis = np.arange(TEST_DPI_MIN, TEST_DPI_MAX, TEST_DPI_STEP)

        # Compute probability for each DPI and DPR
        prob_max = []
        dpi_max = []
        for dpr in test_dprs:
            probs = []
            for dpi in test_dpis:
                prob = 1.0
                for viscode in samples:
                    offset = samples[viscode]['scantime'] * dpr
                    prob *= self.model.get_probability_value(samples[viscode], dpi + offset)
                # prob /= len(samples)
                probs.append(prob)

            arg_max = np.argmax(probs)
            prob_max.append(probs[arg_max])
            dpi_max.append(test_dpis[arg_max])

        # Sanity check
        if np.sum(prob_max) == 0.0:
            print log.WARNING, 'All probabilities equal zero!'
            return None, None

        # Find the DPI and DPR with the highest probability
        arg_max = np.argmax(prob_max)
        return dpi_max[arg_max], test_dprs[arg_max]
