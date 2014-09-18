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
from vgam.progressionmodel import MultiBiomarkerProgressionModel
from vgam.modelfitter import ModelFitter
import fitting.vgam_evaluate_synth as ve


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number_of_test_samples', type=int, default=100, help='the number of test samples')
    parser.add_argument('--recompute_errors', action='store_true', help='recompute the errors of the models')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    parser.add_argument('--number_of_training_samples', type=int, default=1000, help='the number of training samples')
    parser.add_argument('--output_file', type=str, default=None, help='filename of the output image with the plot')
    args = parser.parse_args()

    # Initialise data handler
    data_handler = SynthDataHandler()
    biomarkers = SynthModel.get_biomarker_names()
    viscode_sets = [[0], [0, 1, 2], [0, 1, 2, 3, 4]]

    # Get errors
    errors = get_errors(args, data_handler, biomarkers, viscode_sets)

    # Analyse and plot the test results
    plot_errors(args, errors, biomarkers, viscode_sets)


def get_errors(args, data_handler, biomarkers, viscode_sets):
    evaluation_file = os.path.join(adni.eval_folder, 'eval_synth_fitting_e.p')
    if os.path.isfile(evaluation_file) and not args.recompute_errors:
        # Read test results from file
        print log.INFO, 'Reading test results from {0}...'.format(evaluation_file)
        errors = pickle.load(open(evaluation_file, 'rb'))
    else:
        # Generate test results
        print log.INFO, 'Starting tests...'

        # Initialise model fitter
        model = MultiBiomarkerProgressionModel()
        for biomarker in biomarkers:
            model_file = ve.generate_model(args, data_handler, biomarker, num_samples=args.number_of_training_samples)
            model.add_model(biomarker, model_file)
        fitter = ModelFitter(model)

        # Generate test data and test different biomarker constellations
        num_visits = max([max(vs) for vs in viscode_sets]) + 1
        test_data = ve.generate_test_data(biomarkers, args.number_of_test_samples, num_visits)

        # Compute the errors
        errors = []
        for viscodes in viscode_sets:
            for biomarker in biomarkers:
                errors.append(ve.test_fitting(fitter, test_data, [biomarker], viscodes))
            errors.append(ve.test_fitting(fitter, test_data, biomarkers, viscodes))

        # Save results
        print log.INFO, 'Saving test results to {0}...'.format(evaluation_file)
        pickle.dump(errors, open(evaluation_file, 'wb'))

    return errors


def plot_errors(args, errors, biomarkers, viscode_sets):
    print log.INFO, 'Plotting results...'
    num_viscode_sets = len(viscode_sets)
    num_tests = len(biomarkers) + 1

    # Define figure and grid
    fig, ax1 = plt.subplots(figsize=(10, 6))
    top = 3000
    bottom = -50
    ax1.set_ylim(bottom, top)
    ax1.yaxis.grid(True, linestyle=':', which='major', color='k', alpha=0.4)
    for i in xrange(1, len(viscode_sets)):
        ax1.axvline(i * num_tests + 0.5, linestyle='-', color='k', alpha=0.4)

    # Set labels
    biomarker_strings = {'synth_brain': '$s^{brain}$', 'synth_hipp': '$s^{hipp}$', 'synth_mmse': '$s^{MMSE}$', 'synth_cdrsb': '$s^{CDRSB}$'}
    ax1.set_axisbelow(True)
    ax1.set_title('DPI estimation using different biomarker settings')
    ax1.set_ylabel('Fitting error')
    ax1.set_xticklabels(([biomarker_strings[b] for b in biomarkers] + ['All']) * num_viscode_sets)
    for i in xrange(len(viscode_sets)):
        ax1.text(i * num_tests + 0.5 * (num_tests + 1), -300,
                 '{0} timepoint{1}'.format(len(viscode_sets[i]), '' if len(viscode_sets[i]) == 1 else 's'),
                 horizontalalignment='center')

    # Write mean errors as text
    pos = np.arange(len(errors)) + 1
    medians = [np.median(e) for e in errors]
    upperLabels = [str(np.round(s, 2)) for s in medians]
    weights = (['semibold'] * len(biomarkers) + ['bold']) * num_viscode_sets
    for tick in range(len(errors)):
        ax1.text(pos[tick], 0.95 * top, upperLabels[tick],
                 horizontalalignment='center', size='x-small', weight=weights[tick])

    # Plot boxplots
    boxplot = plt.boxplot(errors, patch_artist=True)
    for i in range(len(boxplot['boxes'])):
        if (i + 1) % num_tests == 0:
            ve.set_boxplot_color(boxplot, i, (0, 0, 1))
        else:
            ve.set_boxplot_color(boxplot, i, (0, 0, 0))

    # Show or save plot
    if args.output_file is not None:
        plt.savefig(args.output_file)
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
