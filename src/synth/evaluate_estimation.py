#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import log as log
from common import synth_tools as st
from common import plotting_tools as pt
from common.datahandler import SynthDataHandler
from common.synthmodel import SynthModel
from common.progressionmodel import ProgressionModel
from common.progressionmodel import MultiBiomarkerProgressionModel
from common.modelfitter import ModelFitter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experiment', type=str, choices=['ex1', 'ex2'], default='ex1', help='the experiment to run')
    parser.add_argument('--number_of_runs', type=int, default=100, help='the number of repeated runs')
    parser.add_argument('--recompute_errors', action='store_true', help='recompute the errors of the models')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    parser.add_argument('--recompute_test_data', action='store_true', help='recompute the test samples')
    parser.add_argument('--number_of_training_samples', type=int, default=1000, help='the number of training samples')
    parser.add_argument('--number_of_test_samples', type=int, default=100, help='the number of test samples')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output image with the plot')
    args = parser.parse_args()

    # Initialise data handler
    biomarkers = SynthModel.get_biomarker_names()

    # Experiment 1
    if args.experiment == 'ex1':
        errors = get_errors_sampling(args, biomarkers)
        plot_errors_sampling(args, biomarkers, errors, args.number_of_training_samples)

    # Experiment 2
    if args.experiment == 'ex2':
        viscode_sets = [[0], [0, 1, 2], [0, 1, 2, 3, 4]]
        biomarker_sets = [[b] for b in biomarkers] + [biomarkers]
        errors = get_errors_data(args, biomarker_sets, viscode_sets)
        plot_errors_data(args, errors, biomarker_sets, viscode_sets)


def get_errors_sampling(args, biomarkers):
    errors = {}
    for biomarker in biomarkers:
        errors.update({biomarker: {}})
        for sampling in ['longitudinal', 'triangular', 'uniform']:
            errors[biomarker].update({sampling: {}})
            e = evaluate_experiment(args, biomarker, sampling)
            errors[biomarker][sampling].update({args.number_of_training_samples: e})
    return errors


def get_errors_data(args, biomarker_sets, viscode_sets):
    # The biomarker sets are the individual biomarkers plus all biomarkers combined
    errors = {}

    for biomarkers in biomarker_sets:
        biomarkers_string = '_'.join(biomarkers)
        errors.update({biomarkers_string: {}})
        for viscodes in viscode_sets:
            viscodes_string = '_'.join([str(v) for v in viscodes])
            if len(biomarkers) == 1:
                e = evaluate_experiment(args, biomarkers[0], 'longitudinal', viscodes=viscodes)
            else:
                e = evaluate_experiment_multiple_biomarkers(args, biomarkers, 'longitudinal', viscodes=viscodes)
            errors[biomarkers_string].update({viscodes_string: e})
    return errors


def evaluate_experiment(args, biomarker, sampling, viscodes=[0]):
    print log.INFO, 'Evaluating {0} model with {1} samples...'.format(biomarker, args.number_of_training_samples)

    num_visits = max(viscodes) + 1
    errors_experiment = []
    for run in xrange(args.number_of_runs):
        model_file = SynthDataHandler.get_model_file(biomarker, num_samples=args.number_of_training_samples,
                                                     sampling=sampling, run=run)
        error_folder = SynthDataHandler.make_dir(SynthDataHandler.get_eval_folder(), biomarker)
        error_file = os.path.join(error_folder, os.path.basename(model_file).replace('.csv', '_test.p'))
        if num_visits > 1:
            error_file = error_file.replace('_test.p', 'v{0}_test.p'.format(num_visits))

        if os.path.isfile(error_file) and not args.recompute_errors:
            print log.SKIP, 'Skipping error computation for {0} samples {1}, run {2}'.format(
                args.number_of_training_samples, sampling, run)
            errors_run = pickle.load(open(error_file, 'rb'))
        else:
            # Generate model
            st.generate_synth_model(biomarker,
                                    recompute_models=args.recompute_models,
                                    num_samples=args.number_of_training_samples,
                                    sampling=sampling, run=run)

            # Initialise fitter
            fitter = ModelFitter(ProgressionModel(biomarker, model_file))

            # Generate test data
            test_data = st.generate_synth_test_data([biomarker],
                                                    args.number_of_test_samples,
                                                    num_visits, run,
                                                    recompute_test_data=args.recompute_test_data)
            errors_run = st.evaluate_synth_fitting(fitter,
                                                   test_data,
                                                   [biomarker],
                                                   viscodes)
            pickle.dump(errors_run, open(error_file, 'wb'))
        errors_experiment.append(np.mean(errors_run))

    return errors_experiment


def evaluate_experiment_multiple_biomarkers(args, biomarkers, sampling, viscodes=[0]):
    print log.INFO, 'Evaluating synthetic models with {0} samples...'.format(args.number_of_training_samples)

    num_visits = max(viscodes) + 1
    errors_experiment = []
    for run in xrange(args.number_of_runs):
        biomarkers_string = '_'.join(biomarkers)
        error_folder = SynthDataHandler.make_dir(SynthDataHandler.get_eval_folder(), biomarkers_string)
        error_file_basename = '{0}_model_{1}_{2}_{3}_test.p'.format(
            biomarkers_string, args.number_of_training_samples, sampling, run)
        error_file = os.path.join(error_folder, error_file_basename)
        if num_visits > 1:
            error_file = error_file.replace('_test.p', 'v{0}_test.p'.format(num_visits))

        if os.path.isfile(error_file) and not args.recompute_errors:
            print log.SKIP, 'Skipping error computation for {0} samples {1}, run {2}'.format(
                args.number_of_training_samples, sampling, run)
            errors_run = pickle.load(open(error_file, 'rb'))
        else:
            # Initialise model fitter
            model = MultiBiomarkerProgressionModel()
            for biomarker in biomarkers:
                model_file = st.generate_synth_model(biomarker,
                                                     recompute_models=args.recompute_models,
                                                     num_samples=args.number_of_training_samples,
                                                     run=run)
                model.add_model(biomarker, model_file)
            fitter = ModelFitter(model)

            # Generate test data
            test_data = st.generate_synth_test_data(biomarkers,
                                                    args.number_of_test_samples,
                                                    num_visits,
                                                    run,
                                                    recompute_test_data=args.recompute_test_data)
            errors_run = st.evaluate_synth_fitting(fitter,
                                                   test_data,
                                                   biomarkers,
                                                   viscodes)
            pickle.dump(errors_run, open(error_file, 'wb'))
        errors_experiment.append(np.mean(errors_run))

    return errors_experiment


def plot_errors_sampling(args, biomarkers, errors, num_samples):
    print log.INFO, 'Plotting error bars...'

    fig, ax = plt.subplots(figsize=(8, 5))
    pt.setup_axes(plt, ax, xgrid=False)

    ax.set_title('Influence of the sampling strategy')
    ax.set_ylabel('Mean progress estimation error')
    ax.set_xticklabels([])

    samplings = ['longitudinal', 'triangular', 'uniform']
    biomarker_strings = {'synth_hipp': '$\mathcal{M}^{HV}$',
                         'synth_mmse': '$\mathcal{M}^{MMSE}$',
                         'synth_cdrsb': '$\mathcal{M}^{CDR-SB}$'}

    # Collect data
    data = []
    medians = []
    for biomarker in biomarkers:
        for sampling in samplings:
            data.append(errors[biomarker][sampling][num_samples])
            medians.append(np.median(errors[biomarker][sampling][num_samples]))

    # Set limits
    max_data = np.max(data)
    ax.set_ylim(0, 1.15 * max_data)

    # Draw boxplot
    boxplot = plt.boxplot(data, patch_artist=True)

    # Set boxplot colours
    colours = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.5), (0.0, 0.5, 0.0)] * len(biomarkers)
    for i in range(len(data)):
        pt.set_boxplot_color(boxplot, i, colours[i])

    # Write median errors as text
    upper_labels = [str(np.round(m, 3)) for m in medians]
    for i in np.arange(len(upper_labels)):
        ax.text(i + 1, 1.02 * max_data, upper_labels[i],
                horizontalalignment='center', size=10, color=colours[i])

    # Write category labels
    for i, biomarker in enumerate(biomarkers):
        plt.text((i + 0.5) * len(samplings) + 0.5, -74,
                 biomarker_strings[biomarker],
                 horizontalalignment='center')

    # Plot legend
    legend = ax.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.0, 0.2), ec=(0.0, 0.0, 0.0, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.5, 0.2), ec=(0.0, 0.0, 0.5, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.5, 0.0, 0.2), ec=(0.0, 0.5, 0.0, 1.0), linewidth=1)],
                       ['{0} sampling'.format(s) for s in samplings], fontsize=10, ncol=3, loc='upper center', framealpha=0.9)
    legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    # Draw horizontal lines
    for x in range(3, len(data), 3):
        plt.axvline(x + 0.5, color='k', alpha=0.4)

    # Show or save plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def plot_errors_data(args, errors, biomarker_sets, viscode_sets):
    print log.INFO, 'Plotting results...'
    num_viscode_sets = len(viscode_sets)
    num_tests = len(biomarker_sets)

    # Define figure and grid
    fig, ax = plt.subplots(figsize=(8, 5))
    pt.setup_axes(plt, ax, xgrid=False)

    ax.set_title('DPI estimation using different biomarker settings')
    ax.set_ylabel('Mean progress estimation error')
    ax.set_xticklabels([])
    biomarker_strings = {'synth_hipp': '$\mathcal{M}^{HV}$',
                         'synth_mmse': '$\mathcal{M}^{MMSE}$',
                         'synth_cdrsb': '$\mathcal{M}^{CDR-SB}$'}

    # Collect data
    data = []
    medians = []
    for viscodes in viscode_sets:
        viscodes_string = '_'.join([str(v) for v in viscodes])
        for biomarkers in biomarker_sets:
            biomarkers_string = '_'.join(biomarkers)
            data.append(errors[biomarkers_string][viscodes_string])
            medians.append(np.median(errors[biomarkers_string][viscodes_string]))

    # Set limits
    max_data = np.max(data)
    ax.set_ylim(0, 1.15 * max_data)

    # Draw boxplot
    boxplot = plt.boxplot(data, patch_artist=True)

    # Set boxplot colours
    colours = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.5), (0.0, 0.5, 0.0), (0.5, 0.0, 0.5)] * num_viscode_sets
    for i in range(len(data)):
        pt.set_boxplot_color(boxplot, i, colours[i])

    # Write median errors as text
    upper_labels = [str(np.round(s, 2)) for s in medians]
    for i in range(len(data)):
        ax.text(i + 1, 1.01 * max_data, upper_labels[i],
                horizontalalignment='center', size=10, color=colours[i])

    # Write category labels
    for i in xrange(len(viscode_sets)):
        ax.text(i * num_tests + 0.5 * (num_tests + 1), -70,
                '{0} timepoint{1}'.format(len(viscode_sets[i]), '' if len(viscode_sets[i]) == 1 else 's'),
                horizontalalignment='center')

    # Plot legend
    legend = ax.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.0, 0.2), ec=(0.0, 0.0, 0.0, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.5, 0.2), ec=(0.0, 0.0, 0.5, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.5, 0.0, 0.2), ec=(0.0, 0.5, 0.0, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.5, 0.0, 0.5, 0.2), ec=(0.5, 0.0, 0.5, 1.0), linewidth=1)],
                       [biomarker_strings[b] for b in biomarker_sets[-1]] + ['All'], fontsize=10, ncol=4, loc='upper center', framealpha=0.9)
    legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    # Draw horizontal lines
    for i in xrange(1, len(viscode_sets)):
        ax.axvline(i * num_tests + 0.5, linestyle='-', color='k', alpha=0.4)

    # Show or save plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
