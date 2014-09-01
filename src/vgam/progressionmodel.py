'''
A class to provide progression modelling functionality.

@author:     Alexander Schmidt-Richberg
@copyright:  2014 Imperial College London. All rights reserved.
@contact:    a.schmidt-richberg@imperial.ac.uk
'''
import math
import numpy as np
import matplotlib.mlab as mlab
from common import log as log


class MultiBiomarkerProgressionModel(object):
    '''
    TODO: classdocs
    '''
    models = {}

    ################################################################################
    #
    # add_model()
    #
    ################################################################################
    def add_model(self, biomarker, model_file, extrapolator='exp'):
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
    def get_probability_value(self, values, progression):
        probability = 1.0
        for biomarker in self.models.keys():
            if biomarker in values:
                value = values[biomarker]
                probability *= self.models[biomarker].get_probability_value(value, progression)
        # probability /= len(self.models)

        return probability


class ProgressionModel(object):
    '''
    TODO: classdocs
    '''

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self, biomarker, model_file, extrapolator='exp'):
        '''
        Initialise the progression model with a data file
        '''
        self.biomarker = biomarker
        self.extrapolator = extrapolator
        self.__intitialise_model(model_file)

    ############################################################################
    #
    # __intitialise_model()
    #
    ############################################################################
    def __intitialise_model(self, model_file):
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

        self.all_progressions = sorted(parameters.keys())
        self.sigma = parameters[self.all_progressions[0]]['sigma']
        self.lmbda = parameters[self.all_progressions[0]]['lambda']
        self.yoffset = parameters[self.all_progressions[0]]['yoffset']

        self.all_mus = [parameters[progr]['mu'] for progr in self.all_progressions]

        self.progressions, self.mus = self.__compute_binned_data(self.all_progressions, self.all_mus)

        self.min_progression = np.min(self.progressions)
        self.max_progression = np.max(self.progressions)

    ############################################################################
    #
    # get_value_range()
    #
    ############################################################################
    def get_value_range(self, std_range=1.5):
        ''' Get the value range estimated by the model in a given std range. '''
        min_curve = self.get_quantile_curve(self.progressions, -std_range)
        max_curve = self.get_quantile_curve(self.progressions, std_range)
        return np.min(min_curve), np.max(max_curve)

    ############################################################################
    #
    # get_quantile_curve()
    #
    ############################################################################
    def get_quantile_curve(self, progressions, std):
        curve = []
        for progression in progressions:
            lmbda, mu, sigma, yoffset = self.get_parameters(progression)
            curve.append(self.__yeojohnson_quantile(lmbda, mu, sigma, std) - yoffset)
        return curve

    ############################################################################
    #
    # get_density_distribution()
    #
    ############################################################################
    def get_density_distribution(self, values, progression):
        lmbda, mu, sigma, yoffset = self.get_parameters(progression)
        probs = [self.__yeojohnson_density(yoffset + v, lmbda, mu, sigma) for v in values]
        return probs

    ############################################################################
    #
    # get_probability_value()
    #
    ############################################################################
    def get_probability_value(self, sample, progression):
        lmbda, mu, sigma, yoffset = self.get_parameters(progression)
        value = sample[self.biomarker] if isinstance(sample, dict) else sample
        return self.__yeojohnson_density(yoffset + value, lmbda, mu, sigma)

    ############################################################################
    #
    # get_parameters()
    #
    ############################################################################
    def get_parameters(self, progression):
        return self.lmbda, self.get_mu(progression), self.sigma, self.yoffset

    ############################################################################
    #
    # get_mu()
    #
    ############################################################################
    def get_mu(self, progression):
        if progression in self.progressions:
            return self.mus[self.progressions.index(progression)]

        elif progression < self.progressions[0]:
            delta_p = self.progressions[1] - self.progressions[0]
            delta_m = self.mus[1] - self.mus[0]

            if self.extrapolator == 'lin':
                # Linear:
                return self.mus[1] - delta_m * ((self.progressions[1] - progression) / delta_p)
            elif self.extrapolator == 'sqrt':
                # Square root:
                return self.mus[1] - delta_m * math.sqrt((self.progressions[1] - progression) / delta_p)
            elif self.extrapolator == 'exp':
                # Exponential:
                index = int(math.floor(len(self.mus) / 2))
                range_m = math.fabs(self.mus[index] - self.mus[0])

                if delta_m > 0:
                    # return self.mus[0] - 1 + math.exp(delta_m / delta_p * (progression - self.progressions[0]))
                    return self.mus[0] - range_m + range_m * math.exp(delta_m / (delta_p * range_m) * (progression - self.progressions[0]))
                else:
                    # return self.mus[0] + 1 - math.exp(-delta_m / delta_p * (progression - self.progressions[0]))
                    return self.mus[0] + range_m - range_m * math.exp(-delta_m / (delta_p * range_m) * (progression - self.progressions[0]))
            else:
                print log.ERROR, 'Unknown extrapolator {0}!'.format(self.extrapolator)
                return 0

        elif progression > self.progressions[-1]:
            delta_p = self.progressions[-1] - self.progressions[-2]
            delta_m = self.mus[-1] - self.mus[-2]

            if self.extrapolator == 'lin':
                # Linear:
                return self.mus[-2] + delta_m * ((progression - self.progressions[-2]) / delta_p)
            elif self.extrapolator == 'sqrt':
                # Square root:
                return self.mus[-2] + delta_m * math.sqrt((progression - self.progressions[-2]) / delta_p)
            elif self.extrapolator == 'exp':
                # Exponential:
                index = int(math.ceil(len(self.mus) / 2))
                range_m = math.fabs(self.mus[index] - self.mus[-1])

                if delta_m > 0:
                    # return self.mus[-1] + 1 - math.exp(-delta_m / delta_p * (progression - self.progressions[-1]))
                    return self.mus[-1] + range_m - range_m * math.exp(-delta_m / (delta_p * range_m) * (progression - self.progressions[-1]))
                else:
                    # return self.mus[-1] - 1 + math.exp(delta_m / delta_p * (progression - self.progressions[-1]))
                    return self.mus[-1] - range_m + range_m * math.exp(delta_m / (delta_p * range_m) * (progression - self.progressions[-1]))
            else:
                print log.ERROR, 'Unknown extrapolator!'
                return 0

        else:
            idx_next = 1
            while self.progressions[idx_next] < progression:
                    idx_next += 1
            idx_prev = idx_next - 1

            progression_prev = self.progressions[idx_prev]
            progression_next = self.progressions[idx_next]
            factor = (progression - progression_prev) / (progression_next - progression_prev)

            mu_prev = self.mus[idx_prev]
            mu_next = self.mus[idx_next]
            return (1 - factor) * mu_prev + factor * mu_next

    ############################################################################
    #
    # __compute_binned_data()
    #
    ############################################################################
    def __compute_binned_data(self, data_x, data_y, bin_size=92):
        bins_x = []
        bins_y = []

        for curr_x in range(np.min(data_x), np.max(data_x), bin_size):
            bin_x = []
            bin_y = []
            for x, y in zip(data_x, data_y):
                if x >= curr_x and x < curr_x + bin_size:
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
    def __yeojohnson_density(self, y, lmbda, mu, sigma):
        ''' Return the probability of a value y given lambda, mu and sigma'''
        return (1 / sigma) * \
            self.__std_normal_dist((self.__yeojohnson(y, lmbda) - mu) / sigma) * \
            math.pow(math.fabs(y) + 1, math.copysign(1, y) * (lmbda - 1))

    ############################################################################
    #
    # __yeojohnson()
    #
    ############################################################################
    def __yeojohnson(self, y, lmbda):
        if y < 0:
            if lmbda == 2:
                return -math.log(-y + 1)
            else:
                return -(math.pow(-y + 1, 2 - lmbda) - 1) / (2 - lmbda)
        else:
            if lmbda == 0:
                return np.log(y + 1)
            else:
                return (math.pow(y + 1, lmbda) - 1) / lmbda

    ############################################################################
    #
    # __std_normal_dist()
    #
    ############################################################################
    def __std_normal_dist(self, x):
        return math.exp(-0.5 * x ** 2) / math.sqrt(2 * math.pi)

    ############################################################################
    #
    # __yeojohnson_quantile()
    #
    ############################################################################
    def __yeojohnson_quantile(self, lmbda, mu, sigma, q):
        ''' Return the value of quantile q'''
        return self.__yeojohnson_inverse(mu + sigma * q, lmbda)

    ############################################################################
    #
    # __yeojohnson_inverse()
    #
    ############################################################################
    def __yeojohnson_inverse(self, x, lmbda):
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