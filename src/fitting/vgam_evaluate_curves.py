#! /usr/bin/env python2.7
import argparse
import os.path
import math
import csv
import numpy as np
from common import log as log
from common import adni_tools as adni
from vgam.datahandler import DataHandler
from vgam.progressionmodel import ProgressionModel
from vgam.modelfitter import ModelFitter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--iteration', type=int, default=0, help='the refinement iteration')
    parser.add_argument('-v', '--value_samples', type=int, default=100, help='the number of values samples')
    parser.add_argument('-p', '--progression_samples', type=int, default=10, help='the number of progression samples')
    parser.add_argument('-s', '--std_range', type=float, default=2.0, help='the range of the standard deviation for the interval computation')
    args = parser.parse_args()

    # Collect data for test
    data_handler = DataHandler(iteration=args.iteration)

    # Compute error for each biomarker
    for biomarker in adni.biomarker_names_essential:
        data_handler.get_model_file(biomarker)
        model_file = data_handler.get_model_file(biomarker)
        if not os.path.isfile(model_file):
            print log.ERROR, 'Model file not found: {0}!'.format(model_file)
        else:
            model = ProgressionModel(biomarker, model_file)
            fitter = ModelFitter(model)

            # Determine value and progression interval
            min_value, max_value = model.get_value_range(std_range=args.std_range)
            values = np.linspace(min_value, max_value, args.value_samples)
            progressions = np.linspace(model.min_progression, model.max_progression, args.progression_samples)
            print log.RESULT, 'Evaluating {0} steps in value interval [{1}, {2}]'.format(args.value_samples, min_value, max_value)
            print log.RESULT, 'Evaluating {0} steps in progression interval [{1}, {2}]'.format(args.progression_samples, model.min_progression, model.max_progression)
            value_step = values[1] - values[0]

            # Compute error
            eval_file = model_file.replace('.csv', '_eval.csv')
            writer = csv.writer(open(eval_file, 'wb'), delimiter=',')
            writer.writerow(['progression', 'error'])

            total_error = 0
            for progression in progressions:
                sample_error = 0
                for value in values:
                    prob_value = model.get_probability_value(value, progression)
                    samples = {'bl': {'scantime': 0, biomarker: value}}
                    estimated_dpi = fitter.get_dpi_for_samples(samples)
                    sample_error += prob_value * np.square(progression - estimated_dpi)
                sample_error = math.sqrt(value_step * sample_error / len(values))
                total_error += sample_error

                writer.writerow([progression, sample_error])
                print log.RESULT, 'Error for progression {0}: {1}'.format(progression, sample_error)

            total_error /= len(model.progressions)
            print log.RESULT, 'Total error: {0}'.format(total_error)


if __name__ == '__main__':
    main()
