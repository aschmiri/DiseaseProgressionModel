#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
from common import log as log
from common import adni_tools as adni
from vgam.datahandler import SynthDataHandler
from vgam.synthmodel import SynthModel
from vgam.progressionmodel import ProgressionModel


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--biomarkers_name', nargs='+', default=None, help='name of the biomarkers to be evaluated')
    parser.add_argument('--recompute_errors', action='store_true', help='recompute the errors of the models')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    parser.add_argument('--experiment_range', type=int, nargs=3, default=[100, 2000, 100], help='the range for the number of samples tested')
    parser.add_argument('--number_of_runs', type=int, default=100, help='the number of repeated runs')
    parser.add_argument('--number_of_progression_steps', type=int, default=10, help='the number of progression steps')
    parser.add_argument('--number_of_value_steps', type=int, default=1000, help='the number of value steps')
    parser.add_argument('--progression_range', type=int, default=2000, help='the width of progression range window used for testing')
    parser.add_argument('--output_file', type=str, default=None, help='filename of the output image with the plot')
    args = parser.parse_args()

    data_handler = SynthDataHandler(args)

    # Experiment 1
    errors = get_errors(args, data_handler)
    plot_errorbars(args, data_handler, errors)
    analyse_errors(args, data_handler, errors)

    # Experiment 2
    num_samples = 1000
    errors = get_errors(args, data_handler, experiments=[num_samples])
    plot_boxplots(args, data_handler, errors, num_samples)


def get_errors(args, data_handler, experiments=None):
    if experiments is None:
        experiments = range(args.experiment_range[0],
                            args.experiment_range[1],
                            args.experiment_range[2])

    errors = {}
    for biomarker in data_handler.get_biomarker_set():
        errors.update({biomarker: {}})
        for sampling in ['longitudinal', 'triangular', 'uniform']:
            errors[biomarker].update({sampling: {}})
            for num_samples in experiments:
                e = evaluate_experiment(args, data_handler, biomarker, sampling, num_samples)
                errors[biomarker][sampling].update({num_samples: e})
    return errors


def evaluate_experiment(args, data_handler, biomarker, sampling, num_samples):
    print log.INFO, 'Evaluating synthetic models with {0} samples...'.format(num_samples)

    errors_experiment = []
    for run in xrange(args.number_of_runs):
        model_file = data_handler.get_model_file(biomarker, num_samples=num_samples, sampling=sampling, run=run)
        error_folder = adni.make_dir(adni.eval_folder, biomarker)
        error_file = os.path.join(error_folder, os.path.basename(model_file).replace('.csv', '.p'))

        if os.path.isfile(error_file) and not args.recompute_errors:
            print log.SKIP, 'Skipping error computation for {0} samples {1}, run {2}'.format(num_samples, sampling, run)
            error_experiment = pickle.load(open(error_file, 'rb'))
        else:
            if os.path.isfile(model_file) and not args.recompute_models:
                print log.SKIP, 'Skipping model generation for {0} samples {1}, run {2}'.format(num_samples, sampling, run)
            else:
                generate_model(args, data_handler, biomarker, sampling, num_samples, run)
            error_experiment = evaluate_model(args, model_file, biomarker)
            pickle.dump(error_experiment, open(error_file, 'wb'))
        errors_experiment.append(error_experiment)

    return errors_experiment


def generate_model(args, data_handler, biomarker, sampling, num_samples, run):
    model_file = data_handler.get_model_file(biomarker)
    model_file_experiment = data_handler.get_model_file(biomarker, num_samples=num_samples, sampling=sampling, run=run)
    samples_file = data_handler.get_samples_file(biomarker)
    samples_file_experiment = data_handler.get_samples_file(biomarker, num_samples=num_samples, sampling=sampling, run=run)

    exec_folder = os.path.join(adni.project_folder, 'src', 'fitting')
    while not os.path.isfile(model_file_experiment):
        call('{0}/vgam_generate_synth_data.py -b {1} -n {2} --sampling {3}'.format(exec_folder, biomarker, num_samples, sampling), shell=True)
        call('{0}/vgam_estimate_curves.py synth -b {1}'.format(exec_folder, biomarker), shell=True)

        if os.path.isfile(model_file) and os.path.isfile(samples_file):
            call(['mv', model_file, model_file_experiment])
            call(['mv', samples_file, samples_file_experiment])
        else:
            print log.WARNING, 'Failed to generate model, retrying...'


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


def plot_errorbars(args, data_handler, errors):
    print log.INFO, 'Plotting error bars...'
    fig = plt.figure()

    linestyle = {'longitudinal': '-', 'triangular': '--', 'uniform': ':'}
    color = {'synth_hipp': 'c', 'synth_brain': 'b', 'synth_mmse': 'g', 'synth_cdrsb': 'r'}

    experiments = range(args.experiment_range[0],
                        args.experiment_range[1],
                        args.experiment_range[2])

    for biomarker in data_handler.get_biomarker_set():
        for sampling in ['longitudinal', 'triangular', 'uniform']:
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

    plt.legend()
    plt.show()

    # Show or save plot
    if args.output_file is not None:
        plt.savefig(args.output_file)
    else:
        plt.show()
    plt.close(fig)


def plot_boxplots(args, data_handler, errors, num_samples):
    print log.INFO, 'Plotting error bars...'
    fig = plt.figure()

    data = []
    labels = []
    for biomarker in data_handler.get_biomarker_set():
        for sampling in ['longitudinal', 'triangular', 'uniform']:
            data.append(errors[biomarker][sampling][num_samples])
            labels.append('{0} {1}'.format(biomarker, sampling))

    plt.boxplot(data)
    plt.xticks(np.arange(len(labels)) + 1, labels)
    plt.title('Comparison of different sampling methods for {0} samples'.format(num_samples))
    plt.grid(True, axis='y', linestyle='--', which='major', color='lightgrey', alpha=0.7)
    plt.legend()

    for x in range(3, len(data), 3):
        plt.axvline(x + 0.5, color='lightgrey', alpha=0.7)

    # Show or save plot
    if args.output_file is not None:
        plt.savefig(args.output_file)
    else:
        plt.show()
    plt.close(fig)


def analyse_errors(args, data_handler, errors):
    experiments = range(args.experiment_range[0],
                        args.experiment_range[1],
                        args.experiment_range[2])

    # Compute mean error uniform vs. triangular
    mean_difference_lon = 0.0
    mean_difference_tri = 0.0
    for biomarker in data_handler.get_biomarker_set():
        error_curve_uniform = []
        error_curve_triangu = []
        error_curve_longitu = []
        for experiment in experiments:
            error_curve_uniform.append(errors[biomarker]['uniform'][experiment])
            error_curve_triangu.append(errors[biomarker]['triangular'][experiment])
            error_curve_longitu.append(errors[biomarker]['longitudinal'][experiment])

        mean_difference_lon += np.mean(np.array(error_curve_longitu) -
                                       np.array(error_curve_uniform))
        mean_difference_tri += np.mean(np.array(error_curve_triangu) -
                                       np.array(error_curve_uniform))
    mean_difference_lon /= len(data_handler.get_biomarker_set())
    mean_difference_tri /= len(data_handler.get_biomarker_set())

    print log.RESULT, 'Mean difference uniform vs. triangular: {0}'.format(mean_difference_tri)
    print log.RESULT, 'Mean difference uniform vs. longitudinal: {0}'.format(mean_difference_lon)


if __name__ == '__main__':
    main()
