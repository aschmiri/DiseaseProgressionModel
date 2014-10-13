#! /usr/bin/env python2.7
import argparse
import os.path
import math
import csv
import joblib as jl
import numpy as np
import matplotlib.mlab as mlab
from common import log as log
from vgam.datahandler import DataHandler
from vgam.progressionmodel import ProgressionModel
from vgam.modelfitter import ModelFitter


def main():
    parser = argparse.ArgumentParser()
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('-v', '--value_samples', type=int, default=100, help='the number of values samples')
    parser.add_argument('-p', '--progress_samples', type=int, default=10, help='the number of progress samples')
    parser.add_argument('-q', '--quantiles', type=float, nargs=2, default=[0.01, 0.99], help='the quantiles for the interval computation')
    parser.add_argument('-n', '--nr_threads', type=int, default=4, help='number of threads')
    args = parser.parse_args()

    # Collect data for test
    data_handler = DataHandler.get_data_handler(args)

    # Compute error for each biomarker
    biomarkers = data_handler.get_biomarker_set()
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(evaluate_biomarker)(args, data_handler, biomarker) for biomarker in biomarkers)

    sort_biomarkers(data_handler, biomarkers)


def evaluate_biomarker(args, data_handler, biomarker):
    data_handler.get_model_file(biomarker)
    model_file = data_handler.get_model_file(biomarker)
    eval_file = model_file.replace('.csv', '_eval.csv')

    if os.path.isfile(eval_file):
        print log.SKIP, 'Evaluation file already existing: {0}'.format(eval_file)
    elif not os.path.isfile(model_file):
        print log.ERROR, 'Model file not found: {0}!'.format(model_file)
    else:
        model = ProgressionModel(biomarker, model_file)
        fitter = ModelFitter(model)

        # Determine value and progress interval
        min_value, max_value = model.get_value_range(quantiles=args.quantiles)
        values = np.linspace(min_value, max_value, args.value_samples)
        progresses = np.linspace(model.min_progress, model.max_progress, args.progress_samples)
        print log.RESULT, 'Evaluating {0} steps in value interval [{1}, {2}]'.format(args.value_samples, min_value, max_value)
        print log.RESULT, 'Evaluating {0} steps in progress interval [{1}, {2}]'.format(args.progress_samples, model.min_progress, model.max_progress)
        value_step = values[1] - values[0]

        # Compute error
        writer = csv.writer(open(eval_file, 'wb'), delimiter=',')
        writer.writerow(['progress', 'error'])

        total_error = 0
        for progress in progresses:
            sample_error = 0
            for value in values:
                prob_value = model.get_probability_value(value, progress)
                samples = {'bl': {'scantime': 0, biomarker: value}}
                estimated_dpi = fitter.get_dpi_for_samples(samples)
                sample_error += prob_value * np.square(progress - estimated_dpi)
            sample_error = math.sqrt(value_step * sample_error / len(values))
            total_error += sample_error

            writer.writerow([progress, sample_error])
            print log.RESULT, 'Error for progress {0}: {1}'.format(progress, sample_error)

        total_error /= len(progresses)
        print log.RESULT, 'Total error: {0}'.format(total_error)


def sort_biomarkers(data_handler, biomarkers):
    errors = []
    for biomarker in biomarkers:
        data_handler.get_model_file(biomarker)
        model_file = data_handler.get_model_file(biomarker)
        eval_file = model_file.replace('.csv', '_eval.csv')

        if os.path.isfile(eval_file):
            m = mlab.csv2rec(eval_file)
            errors.append(np.mean(m['error']))
        else:
            log.ERROR, 'Evaluation file not found: {0}'.format(eval_file)
            errors.append(float('inf'))

    print log.RESULT, 'Sorted biomarkers:'
    max_length = np.max([len(biomarker) for biomarker in biomarkers])
    idx = np.argsort(errors)
    for i in idx:
        print log.RESULT, ('{0:>' + str(max_length) + '}  {1}').format(biomarkers[i], errors[i])
    print log.RESULT, 'Mean error: {0}'.format(np.mean(errors))

if __name__ == '__main__':
    main()
