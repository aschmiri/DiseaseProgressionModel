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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experiment', type=str, choices=['ex1', 'ex2', 'ex3', 'ex4'], default='ex1', help='the experiment to run')
    parser.add_argument('-m', '--metric', type=str, choices=['area', 'peakdist', 'maxdist'], default='area', help='the metric to be evaluated')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarkers to be evaluated')
    parser.add_argument('--recompute_errors', action='store_true', help='recompute the errors of the models')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    parser.add_argument('--sample_numbers_range', type=int, nargs=3, default=[100, 2000, 100], help='the range for the number of samples tested')
    parser.add_argument('--progress_linspace', type=int, nargs=3, default=[-2000, 2000, 100], help='the width of progress range window used for testing')
    parser.add_argument('--number_of_runs', type=int, default=100, help='the number of repeated runs')
    parser.add_argument('--number_of_value_steps', type=int, default=100, help='the number of value steps')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output image with the plot')
    args = parser.parse_args()

    data_handler = SynthDataHandler(biomarkers=args.biomarkers)
    biomarkers = data_handler.get_biomarker_names()

    # Experiment 1
    if args.experiment == 'ex1':
        sample_numbers = range(args.sample_numbers_range[0],
                               args.sample_numbers_range[1],
                               args.sample_numbers_range[2])
        errors = get_errors_samplings(args, biomarkers, sample_numbers)
        plot_error_bars(args, biomarkers, errors)
        analyse_errors(args, biomarkers, errors)

    # Experiment 2
    if args.experiment == 'ex2':
        num_samples = 1000
        errors = get_errors_samplings(args, biomarkers, [num_samples])
        plot_boxplots_samplings(args, biomarkers, errors, num_samples)

    # Experiment 3
    if args.experiment == 'ex3':
        rate_sigmas = [0.0, 0.1, 0.2]
        errors = get_errors_noisy_rates(args, biomarkers, rate_sigmas)
        plot_boxplots_noisy_data(args, biomarkers, errors, rate_sigmas, noise_on_rate=True)

    # Experiment 4
    if args.experiment == 'ex4':
        conversion_sigmas = [0.0, 200.0, 400.0]
        errors = get_errors_noisy_conversion(args, biomarkers, conversion_sigmas)

        plot_boxplots_noisy_data(args, biomarkers, errors, conversion_sigmas, noise_on_rate=False)


def get_errors_samplings(args, biomarkers, sample_numbers):
    errors = {}
    for biomarker in biomarkers:
        errors.update({biomarker: {}})
        for sampling in ['longitudinal', 'triangular', 'uniform']:
            errors[biomarker].update({sampling: {}})
            for num_samples in sample_numbers:
                e = run_experiment(args, biomarker, sampling, num_samples=num_samples)
                errors[biomarker][sampling].update({num_samples: e})
    return errors


def get_errors_noisy_rates(args, biomarkers, rate_sigmas):
    errors = {}
    for biomarker in biomarkers:
        errors.update({biomarker: {}})
        for rate_sigma in rate_sigmas:
            e = run_experiment(args, biomarker, 'longitudinal', rate_sigma=rate_sigma)
            errors[biomarker].update({rate_sigma: e})
    return errors


def get_errors_noisy_conversion(args, biomarkers, conversion_sigmas):
    errors = {}
    for biomarker in biomarkers:
        errors.update({biomarker: {}})
        for conversion_sigma in conversion_sigmas:
            e = run_experiment(args, biomarker, 'longitudinal', conversion_sigma=conversion_sigma)
            errors[biomarker].update({conversion_sigma: e})
    return errors


def run_experiment(args, biomarker, sampling, num_samples=1000, rate_sigma=0.0, conversion_sigma=0.0):
    print log.INFO, 'Evaluating {0} model (sigmas {1} and {2}) with {3} training samples...'.format(
        biomarker, rate_sigma, conversion_sigma, num_samples)

    errors_experiment = []
    for run in xrange(args.number_of_runs):
        data_handler = SynthDataHandler()
        model_file = data_handler.get_model_file(biomarker,
                                                 num_samples=num_samples,
                                                 sampling=sampling,
                                                 rate_sigma=rate_sigma,
                                                 conversion_sigma=conversion_sigma,
                                                 run=run)

        error_folder = SynthDataHandler.make_dir(data_handler.get_eval_folder(), biomarker)
        error_file_ending = '_{0}.p'.format(args.metric)
        error_file = os.path.join(error_folder, os.path.basename(model_file).replace('.csv', error_file_ending))

        if os.path.isfile(error_file) and not args.recompute_errors:
            print log.SKIP, 'Skipping error computation for {0} samples {1}, sigmas {2} and {3}, run {4}'.format(
                num_samples, sampling, rate_sigma, conversion_sigma, run)
            error_experiment = pickle.load(open(error_file, 'rb'))
        else:
            st.generate_synth_model(biomarker,
                                    recompute_models=args.recompute_models,
                                    num_samples=num_samples,
                                    sampling=sampling,
                                    rate_sigma=rate_sigma,
                                    conversion_sigma=conversion_sigma,
                                    run=run)
            error_experiment = st.evaluate_synth_model(model_file,
                                                       biomarker,
                                                       args.progress_linspace,
                                                       args.number_of_value_steps,
                                                       metric=args.metric)
            pickle.dump(error_experiment, open(error_file, 'wb'))
        errors_experiment.append(error_experiment)

    return errors_experiment


def plot_error_bars(args, biomarkers, errors):
    print log.INFO, 'Plotting error bars...'

    fig, ax = plt.subplots(figsize=(8, 5))
    pt.setup_axes(plt, ax)
    ax.set_title('Influence of the number of training samples')
    ax.set_xlabel('Number of samples')
    ax.set_ylabel('Mean area between PDFs')

    color = {'longitudinal': (0.0, 0.0, 0.0), 'triangular': (0.0, 0.0, 0.5), 'uniform': (0.0, 0.5, 0.0)}

    sample_numbers = range(args.sample_numbers_range[0],
                           args.sample_numbers_range[1],
                           args.sample_numbers_range[2])

    for biomarker in biomarkers:
        for sampling in ['longitudinal', 'triangular', 'uniform']:
            curve_median = []
            curve_err_1 = []
            curve_err_2 = []
            for experiment in sample_numbers:
                errors_experiment = errors[biomarker][sampling][experiment]

                median = np.median(errors_experiment)
                curve_median.append(median)
                curve_err_1.append(median - np.percentile(errors_experiment, 25))
                curve_err_2.append(np.percentile(errors_experiment, 75) - median)

            plt.errorbar(sample_numbers, curve_median, yerr=[curve_err_1, curve_err_2],
                         linestyle='-', color=color[sampling],
                         label='{0} sampling'.format(sampling))

    legend = plt.legend(fontsize=10, ncol=3, loc='upper center', framealpha=0.9)
    legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    # Show or save plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def plot_boxplots_samplings(args, biomarkers, errors, num_samples):
    print log.INFO, 'Plotting error bars...'

    fig, ax = plt.subplots(figsize=(8, 5))
    pt.setup_axes(plt, ax, xgrid=False)

    samplings = ['longitudinal', 'triangular', 'uniform']
    biomarker_strings = {'synth_hipp': '$\mathcal{M}^{HV}$',
                         'synth_mmse': '$\mathcal{M}^{MMSE}$',
                         'synth_cdrsb': '$\mathcal{M}^{CDR-SB}$'}
    ylabels = {'area': 'Mean area between PDFs',
               'peakdist': 'Distance between peaks',
               'maxdist': 'Distance between progress maxima'}

    ax.set_title('Influence of the sampling strategy')
    ax.set_ylabel(ylabels[args.metric])
    plt.tick_params(labelbottom='off')

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
        plt.text((i + 0.5) * len(samplings) + 0.5, -0.07 * max_data,
                 biomarker_strings[biomarker],
                 horizontalalignment='center', size=15)

    # Draw horizontal lines
    for x in range(3, len(data), 3):
        plt.axvline(x + 0.5, color='k', alpha=0.4)

    # Plot legend
    legend = ax.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.0, 0.2), ec=(0.0, 0.0, 0.0, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.5, 0.2), ec=(0.0, 0.0, 0.5, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.5, 0.0, 0.2), ec=(0.0, 0.5, 0.0, 1.0), linewidth=1)],
                       ['{0} sampling'.format(s) for s in samplings], fontsize=10, ncol=3, loc='upper center', framealpha=0.9)
    legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    # Show or save plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def plot_boxplots_noisy_data(args, biomarkers, errors, sigmas, noise_on_rate=True):
    print log.INFO, 'Plotting error bars...'

    fig, ax = plt.subplots(figsize=(8, 5))
    pt.setup_axes(plt, ax, xgrid=False)

    biomarker_strings = {'synth_hipp': '$\mathcal{M}^{HV}$',
                         'synth_mmse': '$\mathcal{M}^{MMSE}$',
                         'synth_cdrsb': '$\mathcal{M}^{CDR-SB}$'}
    ylabels = {'area': 'Mean area between PDFs',
               'peakdist': 'Distance between peaks',
               'maxdist': 'Distance between progress maxima'}

    if noise_on_rate:
        ax.set_title('Influence of variations in progression rate')
    else:
        ax.set_title('Influence of variations in point of conversion')
    ax.set_ylabel(ylabels[args.metric])
    plt.tick_params(labelbottom='off')

    # Collect data
    data = []
    medians = []
    for biomarker in biomarkers:
        for sigma in sigmas:
            data.append(errors[biomarker][sigma])
            medians.append(np.mean(errors[biomarker][sigma]))

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
        plt.text((i + 0.5) * len(sigmas) + 0.5, -0.07 * max_data,
                 biomarker_strings[biomarker],
                 horizontalalignment='center', size=15)

    # Draw horizontal lines
    for x in range(3, len(data), 3):
        plt.axvline(x + 0.5, color='k', alpha=0.4)

    # Plot legend
    legend = ax.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.0, 0.2), ec=(0.0, 0.0, 0.0, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.5, 0.2), ec=(0.0, 0.0, 0.5, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.5, 0.0, 0.2), ec=(0.0, 0.5, 0.0, 1.0), linewidth=1)],
                       ['$\sigma={0}$'.format(s) for s in sigmas], fontsize=10, ncol=3, loc='upper center', framealpha=0.9)
    legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    # Show or save plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def analyse_errors(args, biomarkers, errors):
    sample_numbers = range(args.sample_numbers_range[0],
                           args.sample_numbers_range[1],
                           args.sample_numbers_range[2])

    # Compute mean error uniform vs. triangular
    mean_difference_lon = 0.0
    mean_difference_tri = 0.0
    for biomarker in biomarkers:
        error_curve_uniform = []
        error_curve_triangu = []
        error_curve_longitu = []
        for num_samples in sample_numbers:
            error_curve_uniform.append(errors[biomarker]['uniform'][num_samples])
            error_curve_triangu.append(errors[biomarker]['triangular'][num_samples])
            error_curve_longitu.append(errors[biomarker]['longitudinal'][num_samples])

        mean_difference_lon += np.mean(np.abs(np.array(error_curve_longitu) -
                                              np.array(error_curve_uniform)))
        mean_difference_tri += np.mean(np.abs(np.array(error_curve_triangu) -
                                              np.array(error_curve_uniform)))
    mean_difference_lon /= len(biomarkers)
    mean_difference_tri /= len(biomarkers)

    print log.RESULT, 'Mean difference uniform vs. triangular: {0}'.format(mean_difference_tri)
    print log.RESULT, 'Mean difference uniform vs. longitudinal: {0}'.format(mean_difference_lon)


if __name__ == '__main__':
    main()
