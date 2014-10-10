"""
A class to provide progression modelling functionality.

:author:     Alexander Schmidt-Richberg
:copyright:  2014 Imperial College London. All rights reserved.
:contact:    a.schmidt-richberg@imperial.ac.uk
"""
import math
import numpy as np
from common import log as log
from common import adni_tools as adni
from vgam.progressionmodel import ProgressionModel
from vgam.progressionmodel import MultiBiomarkerProgressionModel


class ModelFitter(object):
    """
    TODO: classdocs
    """
    TEST_SCALE_MIN = 0.1
    TEST_SCALE_MAX = 3.0
    TEST_SCALE_STEP = 0.1

    TEST_DPI_MIN = -3500
    TEST_DPI_MAX = 2500
    TEST_DPI_STEP = 10

    TEST_DPR_MIN = 0.0
    TEST_DPR_MAX = 3.0
    TEST_DPR_STEP = 0.1

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self, model):
        """
        Initialise the progression model fitter.

        :param model: the model
        :type model: ProgressionModel or MultiBiomarkerProgressionModel
        """
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
        """ Return the optimally scaled measurements.

        :param dict measurements: the measurements

        :return: the optimally scaled measurements
        :rtype: dict
        """
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
        """ Return the optimal scaling value for a subject given a number of samples
        and a set of biomarkers.

        :param dict samples: the test samples

        :return: the estimated DPI and DPR
        """
        test_scalings = np.arange(self.TEST_SCALE_MIN,
                                  self.TEST_SCALE_MAX,
                                  self.TEST_SCALE_STEP)

        # Compute probability for each scaling
        probs = []
        for scaling in test_scalings:
            prob = math.pow(10.0, len(samples) - 1)
            for viscode in samples:
                scaled_progress = samples[viscode]['progress'] * scaling
                prob *= self.model.get_probability_value(samples[viscode], scaled_progress)
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
        """ Return the estimated DPI of a subject given a number of samples and
        a set of biomarkers.

        :param samples: the test samples

        :return: the estimated DPI
        """
        test_dpis = np.arange(self.TEST_DPI_MIN,
                              self.TEST_DPI_MAX,
                              self.TEST_DPI_STEP)

        # Compute probability for each scaling
        probs = []
        for dpi in test_dpis:
            prob = math.pow(10.0, len(samples) - 1)
            for viscode in samples:
                offset = samples[viscode]['scantime']
                prob *= self.model.get_probability_value(samples[viscode], dpi + offset)
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
        """ Return the estimated DPI and DPR of a subject given a number of samples
        and a set of biomarkers.

        :param samples: the test samples

        :return: the estimated DPI and DPR
        """
        test_dprs = np.arange(self.TEST_DPR_MIN,
                              self.TEST_DPR_MAX,
                              self.TEST_DPR_STEP)
        test_dpis = np.arange(self.TEST_DPI_MIN,
                              self.TEST_DPI_MAX,
                              self.TEST_DPI_STEP)

        # Compute probability for each DPI and DPR
        prob_max = []
        dpi_max = []
        for dpr in test_dprs:
            probs = []
            for dpi in test_dpis:
                prob = math.pow(10.0, len(samples) - 1)
                for viscode in samples:
                    offset = samples[viscode]['scantime'] * dpr
                    prob *= self.model.get_probability_value(samples[viscode], dpi + offset)
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
