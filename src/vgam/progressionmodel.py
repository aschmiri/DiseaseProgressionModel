"""
A class to provide progression modelling functionality.

:author:     Alexander Schmidt-Richberg
:copyright:  2014 Imperial College London. All rights reserved.
:contact:    a.schmidt-richberg@imperial.ac.uk
"""
import math
import numpy as np
import matplotlib.mlab as mlab
from common import log as log


class MultiBiomarkerProgressionModel(object):
    """
    TODO: classdocs
    """
    models = {}

    ################################################################################
    #
    # add_model()
    #
    ################################################################################
    def add_model(self, biomarker, model_file):
        model = ProgressionModel(biomarker, model_file)
        self.models.update({biomarker: model})

    ################################################################################
    #
    # remove_model()
    #
    ################################################################################
    def remove_model(self, biomarker):
        self.models.pop(biomarker, None)

    ############################################################################
    #
    # get_probability_value()
    #
    ############################################################################
    def get_probability_value(self, values, progress):
        probability = math.pow(10.0, len(self.models.keys()) - 1)
        for biomarker in self.models.keys():
            if biomarker in values:
                value = values[biomarker]
                probability *= self.models[biomarker].get_probability_value(value, progress)

        return probability

    ############################################################################
    #
    # get_mean_min_progress()
    #
    ############################################################################
    def get_mean_min_progress(self):
        mean_min_progress = 0.0
        for biomarker in self.models.keys():
            mean_min_progress += self.models[biomarker].min_progress
        return mean_min_progress / len(self.models)

    ############################################################################
    #
    # get_mean_max_progress()
    #
    ############################################################################
    def get_mean_max_progress(self):
        mean_max_progress = 0.0
        for biomarker in self.models.keys():
            mean_max_progress += self.models[biomarker].max_progress
        return mean_max_progress / len(self.models)


class ProgressionModel(object):
    """
    TODO: classdocs
    """

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self, biomarker, model_file, extrapolator='exp'):
        """
        Initialise the progression model with a data file
        """
        self.biomarker = biomarker
        self.extrapolator = extrapolator
        self.__initialise_model(model_file)

    ############################################################################
    #
    # __initialise_model()
    #
    ############################################################################
    def __initialise_model(self, model_file):
        r = mlab.csv2rec(model_file)
        r_sorted = np.sort(r, order=r.dtype.names[2])
        progrs = r_sorted[r.dtype.names[2]]

        parameters = {}
        for i, progr in enumerate(progrs):
            if progr not in parameters:
                parameters.update({progr: {}})
                parameters[progr].update({'sigma': np.exp(r_sorted['logsigma'][i])})
                parameters[progr].update({'lambda': r_sorted['lambda'][i]})
                parameters[progr].update({'mu': r_sorted['mu'][i]})
            parameters[progr].update({'yoffset': r_sorted['fitmiscyoffset'][i]})

        self.all_progresses = sorted(parameters.keys())
        self.sigma = parameters[self.all_progresses[0]]['sigma']
        self.lmbda = parameters[self.all_progresses[0]]['lambda']
        self.yoffset = parameters[self.all_progresses[0]]['yoffset']

        self.all_mus = [parameters[progr]['mu'] for progr in self.all_progresses]

        self.progresses, self.mus = self.__compute_binned_data(self.all_progresses, self.all_mus)

        self.min_progress = np.min(self.progresses)
        self.max_progress = np.max(self.progresses)

    ############################################################################
    #
    # get_value_range()
    #
    ############################################################################
    def get_value_range(self, quantiles=[0.1, 0.9]):
        """ Get the value range estimated by the model in a given std range. """
        min_curve = self.get_quantile_curve(self.progresses, quantiles[0])
        max_curve = self.get_quantile_curve(self.progresses, quantiles[1])
        return np.min(min_curve), np.max(max_curve)

    ############################################################################
    #
    # approximate_quantile()
    #
    ############################################################################
    def approximate_quantile(self, progress, value, step=0.001):
        from scipy import stats
        lmbda, mu, sigma, yoffset = self.get_parameters(progress)
        for quantile in np.arange(step, 1.0, step):
            std = stats.norm.ppf(quantile, loc=0.0, scale=1.0)
            qvalue = self.__yeojohnson_quantile(lmbda, mu, sigma, std) - yoffset
            if qvalue > value:
                return quantile
        return 1.0 - step

    ############################################################################
    #
    # get_quantile_curve()
    #
    ############################################################################
    def get_quantile_curve(self, progresses, quantile):
        curve = []
        for progress in progresses:
            curve.append(self.get_value_at_quantile(progress, quantile))
        return curve

    ############################################################################
    #
    # get_value_at_quantile()
    #
    ############################################################################
    def get_value_at_quantile(self, progress, quantile):
        from scipy import stats
        std = stats.norm.ppf(quantile, loc=0.0, scale=1.0)
        lmbda, mu, sigma, yoffset = self.get_parameters(progress)
        value = self.__yeojohnson_quantile(lmbda, mu, sigma, std) - yoffset
        return value

    ############################################################################
    #
    # get_density_distribution()
    #
    ############################################################################
    def get_density_distribution(self, values, progress):
        lmbda, mu, sigma, yoffset = self.get_parameters(progress)
        probs = [self.__yeojohnson_density(yoffset + v, lmbda, mu, sigma) for v in values]
        return probs

    ############################################################################
    #
    # get_probability_value()
    #
    ############################################################################
    def get_probability_value(self, sample, progress):
        lmbda, mu, sigma, yoffset = self.get_parameters(progress)
        value = sample[self.biomarker] if isinstance(sample, dict) else sample
        return self.__yeojohnson_density(yoffset + value, lmbda, mu, sigma)

    ############################################################################
    #
    # get_parameters()
    #
    ############################################################################
    def get_parameters(self, progress):
        return self.lmbda, self.get_mu(progress), self.sigma, self.yoffset

    ############################################################################
    #
    # get_mu()
    #
    ############################################################################
    def get_mu(self, progress):
        if progress in self.progresses:
            return self.mus[self.progresses.index(progress)]

        elif progress < self.progresses[0]:
            delta_p = self.progresses[1] - self.progresses[0]
            delta_m = self.mus[1] - self.mus[0]

            if self.extrapolator == 'lin':
                # Linear:
                return self.mus[1] - delta_m * ((self.progresses[1] - progress) / delta_p)
            elif self.extrapolator == 'sqrt':
                # Square root:
                return self.mus[1] - delta_m * math.sqrt((self.progresses[1] - progress) / delta_p)
            elif self.extrapolator == 'exp':
                # Exponential:
                range_m = math.fabs(self.mus[0] - self.get_mu(0))

                if delta_m > 0:
                    return self.mus[0] - range_m + range_m * math.exp(delta_m / (delta_p * range_m) * (progress - self.progresses[0]))
                else:
                    return self.mus[0] + range_m - range_m * math.exp(-delta_m / (delta_p * range_m) * (progress - self.progresses[0]))

                # sig_m = math.copysign(1, delta_m)
                # return self.mus[0] - sig_m * (range_m - range_m * math.exp(math.fabs(delta_m) / (delta_p * range_m) * (progress - self.progresses[0])))
            else:
                print log.ERROR, 'Unknown extrapolator {0}!'.format(self.extrapolator)
                return 0

        elif progress > self.progresses[-1]:
            delta_p = self.progresses[-1] - self.progresses[-2]
            delta_m = self.mus[-1] - self.mus[-2]

            if self.extrapolator == 'lin':
                # Linear:
                return self.mus[-2] + delta_m * ((progress - self.progresses[-2]) / delta_p)
            elif self.extrapolator == 'sqrt':
                # Square root:
                return self.mus[-2] + delta_m * math.sqrt((progress - self.progresses[-2]) / delta_p)
            elif self.extrapolator == 'exp':
                # Exponential:
                range_m = math.fabs(self.mus[-1] - self.get_mu(0))

                if delta_m > 0:
                    return self.mus[-1] + range_m - range_m * math.exp(-delta_m / (delta_p * range_m) * (progress - self.progresses[-1]))
                else:
                    return self.mus[-1] - range_m + range_m * math.exp(delta_m / (delta_p * range_m) * (progress - self.progresses[-1]))

                # sig_m = math.copysign(1, delta_m)
                # return self.mus[-1] - sig_m * (-range_m + range_m * math.exp(math.fabs(delta_m) / (delta_p * range_m) * (self.progresses[-1] - progress)))
            else:
                print log.ERROR, 'Unknown extrapolator!'
                return 0

        else:
            idx_next = 1
            while self.progresses[idx_next] < progress:
                    idx_next += 1
            idx_prev = idx_next - 1

            progress_prev = self.progresses[idx_prev]
            progress_next = self.progresses[idx_next]
            factor = (progress - progress_prev) / (progress_next - progress_prev)

            mu_prev = self.mus[idx_prev]
            mu_next = self.mus[idx_next]
            return (1 - factor) * mu_prev + factor * mu_next

    ############################################################################
    #
    # __compute_binned_data()
    #
    ############################################################################
    @staticmethod
    def __compute_binned_data(data_x, data_y, bin_size=92):
        bins_x = []
        bins_y = []

        for curr_x in range(np.min(data_x), np.max(data_x), bin_size):
            bin_x = []
            bin_y = []
            for x, y in zip(data_x, data_y):
                if curr_x <= x < curr_x + bin_size:
                    bin_x.append(x)
                    bin_y.append(y)

            if len(bin_x) > 1:
                selected_x = np.array(bin_x)
                selected_y = np.array(bin_y)

                bins_x.append(np.mean(selected_x))
                bins_y.append(np.mean(selected_y))

        return bins_x, bins_y

    ############################################################################
    #
    # __yeojohnson_density()
    #
    ############################################################################
    @staticmethod
    def __yeojohnson_density(y, lmbda, mu, sigma):
        """ Return the probability of a value y given lambda, mu and sigma"""
        return (1 / sigma) * \
            ProgressionModel.__std_normal_dist((ProgressionModel.__yeojohnson(y, lmbda) - mu) / sigma) * \
            math.pow(math.fabs(y) + 1, math.copysign(1, y) * (lmbda - 1))

    ############################################################################
    #
    # __yeojohnson()
    #
    ############################################################################
    @staticmethod
    def __yeojohnson(y, lmbda):
        if y < 0:
            if lmbda == 2:
                return -math.log(-y + 1)
            else:
                return -(math.pow(-y + 1, 2 - lmbda) - 1) / (2 - lmbda)
        else:
            if lmbda == 0:
                return math.log(y + 1)
            else:
                return (math.pow(y + 1, lmbda) - 1) / lmbda

    ############################################################################
    #
    # __std_normal_dist()
    #
    ############################################################################
    @staticmethod
    def __std_normal_dist(x):
        return math.exp(-0.5 * x ** 2) / math.sqrt(2 * math.pi)

    ############################################################################
    #
    # __yeojohnson_quantile()
    #
    ############################################################################
    @staticmethod
    def __yeojohnson_quantile(lmbda, mu, sigma, q):
        """ Return the value of quantile q"""
        return ProgressionModel.__yeojohnson_inverse(mu + sigma * q, lmbda)

    ############################################################################
    #
    # __yeojohnson_inverse()
    #
    ############################################################################
    @staticmethod
    def __yeojohnson_inverse(x, lmbda):
        if x < 0:
            if lmbda == 2:
                return 1 - math.exp(-x)
            else:
                return 1 - math.pow(-(2 - lmbda) * x + 1, 1 / (2 - lmbda))
        else:
            if lmbda == 0:
                return math.exp(x) - 1
            else:
                return math.pow(lmbda * x + 1, 1 / lmbda) - 1
