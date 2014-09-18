#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
from common import log as log
from common import adni_tools as adni
from vgam.datahandler import SynthDataHandler
from vgam.synthmodel import SynthModel
from vgam.progressionmodel import ProgressionModel
from vgam.modelfitter import ModelFitter
import fitting.vgam_evaluate_synth as ve


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number_of_test_samples', type=int, default=100, help='the number of test samples')
    parser.add_argument('--number_of_runs', type=int, default=10, help='the number of repeated runs')
    parser.add_argument('--recompute_errors', action='store_true', help='recompute the errors of the models')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    parser.add_argument('--number_of_training_samples', type=int, default=1000, help='the number of training samples')
    parser.add_argument('--output_file', type=str, default=None, help='filename of the output image with the plot')
    args = parser.parse_args()

    # Initialise data handler
    data_handler = SynthDataHandler()
    biomarkers = SynthModel.get_biomarker_names()

    # Get errors
    errors = get_errors(args, data_handler, biomarkers)

    # Analyse and plot the test results
    plot_boxplots(args, data_handler, errors, args.number_of_training_samples)


def get_errors(args, data_handler, biomarkers):
    errors = {}
    for biomarker in biomarkers:
        errors.update({biomarker: {}})
        for sampling in ['longitudinal', 'triangular', 'uniform']:
            errors[biomarker].update({sampling: {}})
            e = evaluate_experiment(args, data_handler, biomarker, sampling, args.number_of_training_samples)
            errors[biomarker][sampling].update({args.number_of_training_samples: e})
    return errors


def evaluate_experiment(args, data_handler, biomarker, sampling, num_samples):
    print log.INFO, 'Evaluating synthetic models with {0} samples...'.format(num_samples)

    errors_experiment = []
    for run in xrange(args.number_of_runs):
        model_file = data_handler.get_model_file(biomarker, num_samples=num_samples, sampling=sampling, run=run)
        error_folder = adni.make_dir(adni.eval_folder, biomarker)
        error_file = os.path.join(error_folder, os.path.basename(model_file).replace('.csv', '_test.p'))

        if os.path.isfile(error_file) and not args.recompute_errors:
            print log.SKIP, 'Skipping error computation for {0} samples {1}, run {2}'.format(num_samples, sampling, run)
            error_experiment = pickle.load(open(error_file, 'rb'))
        else:
            # Generate model
            ve.generate_model(args, data_handler, biomarker, num_samples=num_samples, sampling=sampling, run=run)

            # Initialise fitter
            fitter = ModelFitter(ProgressionModel(biomarker, model_file))

            # Generate test data
            test_data = ve.generate_test_data([biomarker], args.number_of_test_samples, 1)
            error_experiment = ve.test_fitting(fitter, test_data, [biomarker])
            pickle.dump(error_experiment, open(error_file, 'wb'))
        errors_experiment.append(error_experiment)

    return errors_experiment


def plot_boxplots(args, data_handler, errors, num_samples):
    print log.INFO, 'Plotting error bars...'
    fig = plt.figure()

    biomarkers = data_handler.get_biomarker_set()
    samplings = ['longitudinal', 'triangular', 'uniform']
    biomarker_strings = {'synth_brain': '$s^{brain}$', 'synth_hipp': '$s^{hipp}$', 'synth_mmse': '$s^{MMSE}$', 'synth_cdrsb': '$s^{CDRSB}$'}
    data = []
    labels = []
    for biomarker in biomarkers:
        for sampling in samplings:
            data.append(errors[biomarker][sampling][num_samples])
            labels.append(sampling)

    for i, biomarker in enumerate(biomarkers):
        plt.text((i + 0.5) * len(samplings) + 0.5, -0.03,
                 biomarker_strings[biomarker],
                 horizontalalignment='center')

    boxplot = plt.boxplot(data, patch_artist=True)
    plt.xticks(np.arange(len(labels)) + 1, labels)
    plt.title('Comparison of different sampling methods for {0} samples'.format(num_samples))
    plt.grid(True, axis='y', linestyle='--', which='major', color='k', alpha=0.4)

    for x in range(3, len(data), 3):
        plt.axvline(x + 0.5, color='k', alpha=0.4)
    for i in range(len(data)):
        ve.set_boxplot_color(boxplot, i, (0, 0, 0))

    # Show or save plot
    if args.output_file is not None:
        plt.savefig(args.output_file)
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
