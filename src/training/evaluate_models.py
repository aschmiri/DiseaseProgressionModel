#! /usr/bin/env python2.7
import argparse
import os.path
import math
import csv
import joblib as jl
import numpy as np
import matplotlib.mlab as mlab
from common import log as log
from common.datahandler import DataHandler
from common.progressionmodel import ProgressionModel
from common.modelfitter import ModelFitter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-p', '--phase', default=None, choices=DataHandler.get_phase_choices(), help='the phase for which the model is to be trained')
    parser.add_argument('-n', '--nr_threads', type=int, default=4, help='number of threads')
    parser.add_argument('--recompute_metric', action='store_true', help='recompute the metric')
    parser.add_argument('--value_samples', type=int, default=100, help='the number of values samples')
    parser.add_argument('--progress_samples', type=int, default=50, help='the number of progress samples')
    parser.add_argument('--quantiles', type=float, nargs=2, default=[0.01, 0.99], help='the quantiles for the interval computation')
    parser.add_argument('--metric', type=str, default='cover', help='the metric used for the evaluation')
    args = parser.parse_args()

    # Collect data for test
    data_handler = DataHandler.get_data_handler(method=args.method,
                                                biomarkers=args.biomarkers,
                                                phase=args.phase)

    # Compute error for each biomarker
    biomarkers = data_handler.get_biomarker_names()
    evaluation_function = evaluate_biomarker_cover if args.metric == 'cover' else evaluate_biomarker_disc
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(evaluation_function)(args, data_handler, biomarker) for biomarker in biomarkers)

    sort_biomarkers(args, data_handler, biomarkers)


def evaluate_biomarker_cover(args, data_handler, biomarker):
    model_file = data_handler.get_model_file(biomarker)
    eval_file = model_file.replace('.csv', '_eval_{0}.csv'.format(args.metric))

    if os.path.isfile(eval_file) and not args.recompute_metric:
        print log.SKIP, 'Evaluation file already existing: {0}'.format(eval_file)
    elif not os.path.isfile(model_file):
        print log.ERROR, 'Model file not found: {0}!'.format(model_file)
    else:
        model = ProgressionModel(biomarker, model_file)

        # Determine value and progress interval
        progresses = np.linspace(model.min_progress, model.max_progress, args.progress_samples)
        median_curve = model.get_quantile_curve(progresses, 0.5)
        min_value = np.min(median_curve)
        max_value = np.max(median_curve)

        print log.INFO, 'Evaluating {0} steps in progress interval [{1}, {2}] for values in [{3}, {4}]...'.format(
            args.progress_samples, progresses[0], progresses[-1], min_value, max_value)

        # Compute error
        writer = csv.writer(open(eval_file, 'wb'), delimiter=',')
        writer.writerow(['progress', 'error'])

        # Compute error
        total_error = 0
        for progress in progresses:
            min_q = model.approximate_quantile(progress, min_value)
            max_q = model.approximate_quantile(progress, max_value)
            quantile_range = max_q - min_q
            total_error += quantile_range

            writer.writerow([progress, quantile_range])

        total_error /= len(progresses)
        print log.RESULT, 'Total error {0}: {1}'.format(biomarker, total_error)


def evaluate_biomarker_disc(args, data_handler, biomarker):
    model_file = data_handler.get_model_file(biomarker)
    eval_file = model_file.replace('.csv', '_eval_{0}.csv'.format(args.metric))

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
        print log.INFO, 'Evaluating {0} steps in value interval [{1}, {2}]...'.format(args.value_samples, min_value, max_value)
        print log.INFO, 'Evaluating {0} steps in progress interval [{1}, {2}]...'.format(args.progress_samples, model.min_progress, model.max_progress)
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
                estimated_dpi = fitter.get_dpi_for_samples(samples, phase=args.args)
                sample_error += prob_value * np.square(progress - estimated_dpi)
            sample_error = math.sqrt(value_step * sample_error / len(values))
            total_error += sample_error

            writer.writerow([progress, sample_error])
            print log.RESULT, 'Error for progress {0}: {1}'.format(progress, sample_error)

        total_error /= len(progresses)
        print log.RESULT, 'Total error: {0}'.format(total_error)


def sort_biomarkers(args, data_handler, biomarkers):
    errors = []
    for biomarker in biomarkers:
        model_file = data_handler.get_model_file(biomarker)
        eval_file = model_file.replace('.csv', '_eval_{0}.csv'.format(args.metric))

        if os.path.isfile(eval_file):
            m = mlab.csv2rec(eval_file)
            errors.append(np.mean(m['error']))
        else:
            log.ERROR, 'Evaluation file not found: {0}'.format(eval_file)
            errors.append(float(0.0))

    print log.RESULT, 'Sorted biomarkers:'
    max_length = np.max([len(biomarker) for biomarker in biomarkers])
    idx = np.argsort(errors)
    for i in reversed(idx):
        print log.RESULT, ('{0:>' + str(max_length) + '}  {1}').format(biomarkers[i], errors[i])
    print log.RESULT, 'Mean error: {0}'.format(np.mean(errors))

if __name__ == '__main__':
    main()
