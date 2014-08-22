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
from common import adni_tools as adni
from common import adni_plot as aplt
from vgam.progressionmodel import ProgressionModel
from vgam.datahandler import DataHandler


def main():
    parser = argparse.ArgumentParser()
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('-e', '--extrapolator', type=str, choices=['lin', 'sqrt', 'exp'], default='exp', help='the type of extrapolator')
    parser.add_argument('-p', '--no_points', action='store_true', default=False, help='do not plot points')
    parser.add_argument('-d', '--no_densities', action='store_true', default=False, help='do not plot densities')
    parser.add_argument('--plot_mu', action='store_true', default=False, help='plot mu')
    parser.add_argument('--plot_errors', action='store_true', default=False, help='plot th errors')
    parser.add_argument('--save_file', action='store_true', default=False, help='save the plots as a file')
    args = parser.parse_args()

    data_handler = DataHandler(args)
    for biomarker in data_handler.get_biomarker_set():
        plot_model(args, data_handler, biomarker)


def plot_model(args, data_handler, biomarker):
    samples_file = data_handler.get_samples_file(biomarker)
    model_file = data_handler.get_model_file(biomarker)
    if not os.path.isfile(samples_file) or not os.path.isfile(model_file):
        print log.ERROR, 'Files not available for {0}!'.format(biomarker)
        return

    print log.INFO, 'Generating plot for {0}...'.format(biomarker)

    #
    # Read model
    #
    pm = ProgressionModel(biomarker, model_file, extrapolator=args.extrapolator)
    min_progression_extra = pm.min_progression - 0.3 * (pm.max_progression - pm.min_progression)
    max_progression_extra = pm.max_progression + 0.3 * (pm.max_progression - pm.min_progression)
    progression_linspace = np.linspace(min_progression_extra, max_progression_extra, 100)

    #
    # Setup plot
    #
    plt.figure(figsize=(12, 5), dpi=100)
    if not args.no_densities:
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
    else:
        ax1 = plt.subplot(1, 1, 1)

    ax1.set_title('Percentile curves for {0}'.format(biomarker))
    ax1.set_xlabel('Disease progression relative to point of conversion')
    ax1.set_ylabel('Volume' if biomarker in adni.volume_names else 'Score')
    ax1.set_xlim(min_progression_extra, max_progression_extra)
    ax1.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.8, 0.8, 0.0), linewidth=0),
                mpl.patches.Rectangle((0, 0), 1, 1, fc=(1.0, 0.0, 0.0), linewidth=0)],
               ['MCI', 'AD'], fontsize=10)

    #
    # Plot the percentile curves
    #
    ax1.axvline(pm.min_progression, color='0.15', linestyle=':', alpha=0.8)
    ax1.axvline(pm.max_progression, color='0.15', linestyle=':', alpha=0.8)

    stds = [-1.5, -1.0, -0.5, 0, +0.5, +1.0, +1.5]
    greyvals = ['0.45', '0.3', '0.15', '0', '0.15', '0.3', '0.45']
    min_vals = []
    max_vals = []
    for (greyval, std) in zip(greyvals, stds):
        curve = pm.get_quantile_curve(progression_linspace, std)
        min_vals.append(np.min(curve))
        max_vals.append(np.max(curve))
        ax1.plot(progression_linspace, curve, color=greyval)

        label = '${0} \sigma$'.format(std)
        ax1.text(pm.progressions[-1], curve[-1], label, fontsize=11)

    min_val = np.min(min_vals)
    max_val = np.max(max_vals)

    #
    # Plot parameter mu
    #
    if args.plot_mu:
        # Get second axis of plot 1
        ax1b = ax1.twinx()

        # Plot all progressions
        ax1b.scatter(pm.all_progressions, pm.all_mus, color='b', marker='o', linewidths=0, alpha=0.2)
        ax1b.text(pm.progressions[-1], pm.mus[-1], '$\mu$', color='b', fontsize=11)

        # Plot binned progressions
        ax1b.scatter(pm.progressions, pm.mus, color='b', marker='x')

        # Plot interpolated model
        mus = [pm.get_mu(p) for p in progression_linspace]
        ax1b.plot(progression_linspace, mus, color='b')
        ax1b.set_xlim(min_progression_extra, max_progression_extra)

    #
    # Plot errors
    #
    if args.plot_errors:
        eval_file = model_file.replace('.csv', '_eval.csv')
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
    if not args.no_points:
        m = mlab.csv2rec(samples_file)
        progr_points = m['progress']
        value_points = m['value']
        diagn_points = [0.5 if p < 0 else 1.0 for p in progr_points]

        print log.INFO, 'Plotting {0} sample points...'.format(len(progr_points))
        ax1.scatter(progr_points, value_points, c=diagn_points, vmin=0.0, vmax=1.0, linewidths=0, cmap=aplt.progression_cmap, alpha=0.25)

    #
    # Plot PDFs
    #
    if not args.no_densities:
        ax2.set_title('Probability density function for {0}'.format(biomarker))
        ax2.set_xlabel('Volume' if biomarker in adni.volume_names else 'Score')
        ax2.set_ylabel('Probability')

        values = np.linspace(min_val, max_val, 250)
        progr_samples = [-2000, -1500, -1000, -500, 0, 500, 1000, 1500, 2000]

        sample_cmap = cmx.ScalarMappable(
            norm=colors.Normalize(vmin=-len(progr_samples) + 1, vmax=(len(progr_samples) - 1)),
            cmap=plt.get_cmap(aplt.progression_cmap))

        for progr in progr_samples:
            sample_color = sample_cmap.to_rgba(progr_samples.index(progr))
            ax1.axvline(progr, color=sample_color, linestyle='--', alpha=0.3)
            ax2.set_xlim(min_val, max_val)

            linestyle = '--' if progr < pm.min_progression or progr > pm.max_progression else '-'
            probs = pm.get_density_distribution(values, progr)
            ax2.plot(values, probs, label=str(progr), color=sample_color, linestyle=linestyle)

        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(handles2, labels2, fontsize=10)

    #
    # Draw or save the plot
    #
    plt.tight_layout()
    if args.save_file:
        plot_filename = model_file.replace('.csv', '.pdf')
        plt.savefig(plot_filename, dpi=100)
    else:
        plt.show()

if __name__ == '__main__':
    main()
