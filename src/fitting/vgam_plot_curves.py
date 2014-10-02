#! /usr/bin/env python2.7
import os.path
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cmx
from common import log as log
from common import adni_plot as aplt
from vgam.progressionmodel import ProgressionModel
from vgam.datahandler import DataHandler
from vgam.synthmodel import SynthModel
import fitting.vgam_evaluation as ve


def main():
    parser = argparse.ArgumentParser()
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('-e', '--extrapolator', type=str, choices=['lin', 'sqrt', 'exp'], default='exp', help='the type of extrapolator')
    parser.add_argument('--no_model', action='store_true', default=False, help='do not plot the fitted model')
    parser.add_argument('--no_points', action='store_true', default=False, help='do not plot points')
    parser.add_argument('--points_alpha', type=float, default=0.25, help='alpha value of the plotted points')
    parser.add_argument('--no_densities', action='store_true', default=False, help='do not plot densities')
    parser.add_argument('--no_sample_lines', action='store_true', default=False, help='do not plot the sample lines')
    parser.add_argument('--only_densities', action='store_true', default=False, help='only plot densities')
    parser.add_argument('--no_extrapolation', action='store_true', default=False, help='do not extrapolate the model')
    parser.add_argument('--plot_mu', action='store_true', default=False, help='plot mu')
    parser.add_argument('--plot_errors', action='store_true', default=False, help='plot the errors')
    parser.add_argument('--plot_synth_model', action='store_true', default=False, help='plot density distributions for synthetic data')
    parser.add_argument('--save_file', action='store_true', default=False, help='save the plots as a file')
    parser.add_argument('--output_file', type=str, default=None, help='filename of the output file')
    args = parser.parse_args()

    data_handler = DataHandler.get_data_handler(args)
    for biomarker in data_handler.get_biomarker_set():
        plot_model(args, data_handler, biomarker)


def plot_model(args, data_handler, biomarker):
    model_file = data_handler.get_model_file(biomarker)
    if not os.path.isfile(model_file):
        print log.ERROR, 'Model file not found: {0}'.format(model_file)
        return

    print log.INFO, 'Generating plot for {0}...'.format(biomarker)
    plot_synth_model = args.plot_synth_model and biomarker in SynthModel.get_biomarker_names()

    #
    # Read model
    #
    pm = ProgressionModel(biomarker, model_file, extrapolator=args.extrapolator)
    progression_extrapolate = 0.3 * (pm.max_progression - pm.min_progression)
    min_progression_extrapolate = pm.min_progression - progression_extrapolate
    max_progression_extrapolate = pm.max_progression + progression_extrapolate
    progression_linspace_ex1 = np.linspace(min_progression_extrapolate, pm.min_progression, 20)
    progression_linspace_int = np.linspace(pm.min_progression, pm.max_progression, 60)
    progression_linspace_ex2 = np.linspace(pm.max_progression, max_progression_extrapolate, 20)

    # Calc min and max val in interval [-1.5 sigma, +1.5 sigma]
    progression_linspace = np.linspace(min_progression_extrapolate, max_progression_extrapolate, 100)
    min_val = float('inf')
    max_val = float('-inf')
    for std in [-1.0, 1.0]:
        curve = pm.get_quantile_curve(progression_linspace, std)
        min_val = min(min_val, np.min(curve))
        max_val = max(max_val, np.max(curve))

    #
    # Setup plot
    #
    figure_width = 6 if args.no_densities or args.only_densities else 12
    fig = plt.figure(figsize=(figure_width, 5))
    if args.only_densities:
        ax1 = None
        ax2 = plt.subplot(1, 1, 1)
        ve.setup_axes(plt, ax2)
    elif args.no_densities:
        ax1 = plt.subplot(1, 1, 1)
        ax2 = None
        ve.setup_axes(plt, ax1)
    else:
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
        ve.setup_axes(plt, ax1)
        ve.setup_axes(plt, ax2)

    if not args.only_densities:
        ax1.set_title('Percentile curves for {0}'.format(biomarker))
        ax1.set_xlabel('Disease progression relative to point of conversion')
        ax1.set_ylabel(ve.get_metric_unit(biomarker))
        ax1.set_xlim(min_progression_extrapolate, max_progression_extrapolate)

    #
    # Plot the percentile curves of the fitted model
    #
    if not args.no_model and not args.only_densities:
        ax1.axvline(pm.min_progression, color='0.15', linestyle=':', alpha=0.8)
        ax1.axvline(pm.max_progression, color='0.15', linestyle=':', alpha=0.8)

        stds = [-1.5, -1.0, -0.5, 0, +0.5, +1.0, +1.5]
        grey_values = ['0.6', '0.4', '0.2', '0', '0.2', '0.4', '0.6']
        for grey_value, std in zip(grey_values, stds):
            curve_int = pm.get_quantile_curve(progression_linspace_int, std)
            ax1.plot(progression_linspace_int, curve_int, color=grey_value)

            if not args.no_extrapolation:
                curve_ex1 = pm.get_quantile_curve(progression_linspace_ex1, std)
                curve_ex2 = pm.get_quantile_curve(progression_linspace_ex2, std)
                ax1.plot(progression_linspace_ex1, curve_ex1, '--', color=grey_value)
                ax1.plot(progression_linspace_ex2, curve_ex2, '--', color=grey_value)

            # TODO: label = '${0} \sigma$'.format(std)
            # ax1.text(pm.progressions[-1], curve_2[-1], label, fontsize=11)

    #
    # Plot synthetic model curve
    #
    if plot_synth_model:
        progression_linspace_synth = np.linspace(-2500, 2500, 100)
        quantiles = [0.1, 0.25, 0.5, 0.75, 0.9]
        alphas = [0.4, 0.7, 1.0, 0.7, 0.4]
        for quantile, alpha in zip(quantiles, alphas):
            curve_synth = [SynthModel.get_distributed_value(biomarker, p, cdf=quantile) for p in progression_linspace_synth]
            ax1.plot(progression_linspace_synth, curve_synth, color='b', alpha=alpha)

    #
    # Plot parameter mu
    #
    if args.plot_mu and not args.only_densities:
        # Get second axis of plot 1
        ax1b = ax1.twinx()

        # Plot all progressions
        ax1b.scatter(pm.all_progressions, pm.all_mus, facecolor='b', marker='o', edgecolor='none', alpha=0.2)
        ax1b.text(pm.progressions[-1], pm.mus[-1], '$\mu$', color='b', fontsize=11)

        # Plot binned progressions
        ax1b.scatter(pm.progressions, pm.mus, color='b', marker='x')

        # Plot interpolated model
        mus = [pm.get_mu(p) for p in progression_linspace_int]
        ax1b.plot(progression_linspace_int, mus, color='b')

        if not args.no_extrapolation:
            mus = [pm.get_mu(p) for p in progression_linspace_ex1]
            ax1b.plot(progression_linspace_ex1, mus, '--', color='b')
            mus = [pm.get_mu(p) for p in progression_linspace_ex2]
            ax1b.plot(progression_linspace_ex2, mus, '--', color='b')

        ax1b.set_xlim(min_progression_extrapolate, max_progression_extrapolate)

    #
    # Plot errors
    #
    if args.plot_errors and not args.only_densities:
        eval_file = model_file.replace('.csv', '_eval.csv')
        if not os.path.isfile(eval_file):
            print log.ERROR, 'Evaluation file not found: {0}'.format(eval_file)
        else:
            m = mlab.csv2rec(eval_file)
            progressions = m['progression']
            errors = m['error']

            # Get second axis of plot 1
            ax1b = ax1.twinx()
            ax1b.set_ylim(0, max(150, 1.2 * np.max(errors)))
            ax1b.plot(progressions, errors, color='g', marker='x')
            ax1b.text(progressions[-1], errors[-1], 'Discr.', color='g', fontsize=11)
            ax1b.axhline(np.mean(errors), color='g', linestyle='--', alpha=0.5)

    #
    # Plot points
    #
    if not args.no_points and not args.only_densities:
        samples_file = data_handler.get_samples_file(biomarker)
        if not os.path.isfile(samples_file):
            print log.ERROR, 'Samples file not found: {0}'.format(samples_file)
        else:
            m = mlab.csv2rec(samples_file)
            progr_points = m['progress']
            value_points = m['value']
            diagn_points = [0.5 if p < 0 else 1.0 for p in progr_points]

            print log.INFO, 'Plotting {0} sample points...'.format(len(progr_points))
            ax1.scatter(progr_points, value_points, c=diagn_points, edgecolor='none',
                        vmin=0.0, vmax=1.0, cmap=aplt.progression_cmap, alpha=args.points_alpha)
            ax1.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.8, 0.8, 0.0), linewidth=0),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(1.0, 0.0, 0.0), linewidth=0)],
                       ['MCI', 'AD'], fontsize=10)

    #
    # Plot PDFs
    #
    progr_samples = [-2000, -1500, -1000, -500, 0, 500, 1000, 1500, 2000]
    sample_cmap = cmx.ScalarMappable(
        norm=colors.Normalize(vmin=-len(progr_samples) + 1, vmax=(len(progr_samples) - 1)),
        cmap=plt.get_cmap(aplt.progression_cmap))

    if not args.no_sample_lines and not args.only_densities:
        for progr in progr_samples:
            if not args.no_extrapolation or pm.min_progression < progr < pm.max_progression:
                sample_color = sample_cmap.to_rgba(progr_samples.index(progr))
                linestyle = '--' if progr < pm.min_progression or progr > pm.max_progression else '-'
                ax1.axvline(progr, color=sample_color, linestyle=linestyle, alpha=0.3)

    if not args.no_densities:
        ax2.set_title('Probability density function for {0}'.format(biomarker))
        ax2.set_xlabel(ve.get_metric_unit(biomarker))
        ax2.set_ylabel('Probability')

        values = np.linspace(min_val, max_val, 250)
        for progr in progr_samples:
            if not args.no_extrapolation or pm.min_progression < progr < pm.max_progression:
                sample_color = sample_cmap.to_rgba(progr_samples.index(progr))
                linestyle = '--' if progr < pm.min_progression or progr > pm.max_progression else '-'
                probs = pm.get_density_distribution(values, progr)
                ax2.set_xlim(min_val, max_val)
                ax2.plot(values, probs, label=str(progr), color=sample_color, linestyle=linestyle)

                if plot_synth_model:
                    probs = [SynthModel.get_probability(biomarker, progr, v) for v in values]
                    ax2.plot(values, probs, color='b', linestyle='--')

        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(handles2, labels2, fontsize=10)

    #
    # Draw or save the plot
    #
    plt.tight_layout()
    if args.save_file or args.output_file is not None:
        if args.output_file is not None:
            plot_filename = args.output_file
        else:
            plot_filename = model_file.replace('.csv', '.pdf')
        plt.savefig(plot_filename, transparent=True)
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
