#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
import numpy as np
from subprocess import call
import matplotlib.pyplot as plt
from common import log as log
from common import adni_tools as adni
from vgam.datahandler import SynthDataHandler
from vgam.synthmodel import SynthModel
from vgam.progressionmodel import ProgressionModel


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--biomarker_name', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('--recompute_errors', action='store_true', help='recompute the errors of the models')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    parser.add_argument('--experiment_range', type=int, nargs=3, default=[100, 2000, 100], help='the range for the number of samples tested')
    parser.add_argument('--number_of_runs', type=int, default=100, help='the number of repeated runs')
    parser.add_argument('--number_of_progression_steps', type=int, default=10, help='the number of progression steps')
    parser.add_argument('--number_of_value_steps', type=int, default=1000, help='the number of value steps')
    parser.add_argument('--progression_range', type=int, default=2000, help='the width of progression range window used for testing')
    args = parser.parse_args()

    data_handler = SynthDataHandler(args)

    if args.recompute_errors:
        errors = compute_errors(args, data_handler)
        save_errors(args, errors)
    else:
        errors = read_errors(args)

    plot_errors(args, data_handler, errors)
    analyse_errors(args, data_handler, errors)


def compute_errors(args, data_handler):
    errors = {}
    for biomarker in data_handler.get_biomarker_set():
        errors.update({biomarker: {}})
        for sampling in ['triangular', 'uniform']:
            errors[biomarker].update({sampling: {}})
            for num_samples in xrange(args.experiment_range[0],
                                      args.experiment_range[1],
                                      args.experiment_range[2]):
                e = evaluate_experiment(args, data_handler, biomarker, sampling, num_samples)
                errors[biomarker][sampling].update({num_samples: e})
    return errors


def evaluate_experiment(args, data_handler, biomarker, sampling, num_samples):
    print log.INFO, 'Evaluating synthetic models with {0} samples...'.format(num_samples)

    errors_experiment = []
    for run in xrange(args.number_of_runs):
        model_file = generate_models(args, data_handler, biomarker, sampling, num_samples, run)
        errors_experiment.append(evaluate_model(args, model_file, biomarker))

    return errors_experiment


def generate_models(args, data_handler, biomarker, sampling, num_samples, run):
    experiment_ending = '_{0}_{1}_{2}.csv'.format(num_samples, sampling, run)
    model_file = data_handler.get_model_file(biomarker)
    model_file_experiment = model_file.replace('.csv', experiment_ending)
    samples_file = data_handler.get_samples_file(biomarker)
    samples_file_experiment = samples_file.replace('.csv', experiment_ending)

    if os.path.isfile(model_file_experiment) and os.path.isfile(samples_file_experiment) and not args.recompute_models:
        print log.SKIP, 'Skipping model generation for {0} samples {1}, run {2}'.format(num_samples, sampling, run)
    else:
        exec_folder = os.path.join(adni.project_folder, 'src', 'fitting')
        sampling_arg = '--uniform_progression' if sampling == 'uniform' else ''
        call('{0}/vgam_generate_synth_data.py -b {1} -n {2} {3}'.format(exec_folder, biomarker, num_samples, sampling_arg), shell=True)
        call('{0}/vgam_estimate_curves.py synth -b {1}'.format(exec_folder, biomarker), shell=True)
        call(['mv', model_file, model_file_experiment])
        call(['mv', samples_file, samples_file_experiment])

    return model_file_experiment


def evaluate_model(args, model_file, biomarker):

    # Define progression steps
    pm = ProgressionModel(biomarker, model_file)
    progressions = np.linspace(-args.progression_range,
                               args.progression_range,
                               args.number_of_progression_steps)

    # Define value steps
    min_val = float('inf')
    max_val = float('-inf')
    for std in [-2.0, 2.0]:
        curve = pm.get_quantile_curve(progressions, std)
        min_val = min(min_val, np.min(curve))
        max_val = max(max_val, np.max(curve))
    values = np.linspace(min_val, max_val, args.number_of_value_steps)

    # Get mean error
    error = 0
    for progr in progressions:
        probs_fit = pm.get_density_distribution(values, progr)
        probs_model = [SynthModel.get_probability(biomarker, progr, v) for v in values]
        error += np.sum(np.abs(np.array(probs_fit) - np.array(probs_model)))
    error *= (values[1] - values[0]) / len(progressions)

    return error


def read_errors(args):
    if args.biomarker_name is None:
        file_name = 'eval_synth_model.p'
    else:
        file_name = 'eval_synth_model_{0}.p'.format(args.biomarker_name)
    evaluation_file = os.path.join(adni.eval_folder, file_name)
    return pickle.load(open(evaluation_file, 'rb'))


def save_errors(args, errors):
    if args.biomarker_name is None:
        file_name = 'eval_synth_model.p'
    else:
        file_name = 'eval_synth_model_{0}.p'.format(args.biomarker_name)
    evaluation_file = os.path.join(adni.eval_folder, file_name)
    pickle.dump(errors, open(evaluation_file, 'wb'))


def plot_errors(args, data_handler, errors):
    print log.INFO, 'Plotting error bars...'

    plt.figure()
    linestyle = {'triangular': '-', 'uniform': '--'}
    color = {'synth0': 'c', 'synth1': 'b', 'synth2': 'g', 'synth3': 'r', 'synth4': 'k'}

    experiments = range(args.experiment_range[0],
                        args.experiment_range[1],
                        args.experiment_range[2])

    for biomarker in data_handler.get_biomarker_set():
        for sampling in ['triangular', 'uniform']:
            curve_median = []
            curve_err_1 = []
            curve_err_2 = []
            for experiment in experiments:
                errors_experiment = errors[biomarker][sampling][experiment]

                median = np.median(errors_experiment)
                curve_median.append(median)
                curve_err_1.append(median - np.percentile(errors_experiment, 25))
                curve_err_2.append(np.percentile(errors_experiment, 75) - median)

            plt.errorbar(experiments, curve_median, yerr=[curve_err_1, curve_err_2],
                         linestyle=linestyle[sampling], color=color[biomarker],
                         label='{0} {1}'.format(biomarker, sampling))

    plt.xlim(90, 1010)
    plt.legend()
    plt.show()


def analyse_errors(args, data_handler, errors):
    experiments = range(args.experiment_range[0],
                        args.experiment_range[1],
                        args.experiment_range[2])

    # Compute mean error uniform vs. triangular
    mean_difference = 0.0
    for biomarker in data_handler.get_biomarker_set():
        error_curve_uniform = []
        error_curve_triangu = []
        for experiment in experiments:
            error_curve_uniform.append(errors[biomarker]['uniform'][experiment])
            error_curve_triangu.append(errors[biomarker]['triangular'][experiment])

        mean_difference += np.mean(np.array(error_curve_triangu) -
                                   np.array(error_curve_uniform))
    mean_difference /= len(data_handler.get_biomarker_set())

    print log.RESULT, 'Mean difference uniform vs. triangular: {0}'.format(mean_difference)


if __name__ == '__main__':
    main()
