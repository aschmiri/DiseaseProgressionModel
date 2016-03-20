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
from common import plotting_tools as pt
from common.progressionmodel import ProgressionModel
from common.datahandler import DataHandler
from common.synthmodel import SynthModel


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-p', '--phase', default='mciad', choices=DataHandler.get_phase_choices(), help='the phase for which the model is to be trained')
    parser.add_argument('-e', '--extrapolator', type=str, choices=['lin', 'sqrt', 'exp'], default='exp', help='the type of extrapolator')
    parser.add_argument('--xlim', type=float, nargs=2, default=None, help='force certain x limits for plotting')
    parser.add_argument('--ylim', type=float, nargs=2, default=None, help='force certain y limits for plotting')
    parser.add_argument('--no_model', action='store_true', default=False, help='do not plot the fitted model')
    parser.add_argument('--no_points', action='store_true', default=False, help='do not plot points')
    parser.add_argument('--points_alpha', type=float, default=0.25, help='alpha value of the plotted points')
    parser.add_argument('--no_densities', action='store_true', default=False, help='do not plot densities')
    parser.add_argument('--no_sample_lines', action='store_true', default=False, help='do not plot the sample lines')
    parser.add_argument('--only_densities', action='store_true', default=False, help='only plot densities')
    parser.add_argument('--no_extrapolation', action='store_true', default=False, help='do not extrapolate the model')
    parser.add_argument('--plot_eta', type=str, choices=['lambda', 'mu', 'sigma'], default=None, help='plot a predictor function')
    parser.add_argument('--plot_errors', action='store_true', default=False, help='plot the errors')
    parser.add_argument('--plot_synth_model', action='store_true', default=False, help='plot density distributions for synthetic data')
    parser.add_argument('--plot_quantile_label', action='store_true', default=False, help='plot labels on the quantile curces')
    parser.add_argument('--plot_donohue', action='store_true', default=False, help='plot the trajectory estimated with Donohue et al.')
    parser.add_argument('--save_plots', action='store_true', default=False, help='save the plots with a default filename')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output file')
    args = parser.parse_args()

    data_handler = DataHandler.get_data_handler(method=args.method,
                                                biomarkers=args.biomarkers,
                                                phase=args.phase)
    for biomarker in data_handler.get_biomarker_names():
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
    progress_extrapolate = 0.3 * (pm.max_progress - pm.min_progress)
    min_progress_extrapolate = int(pm.min_progress - progress_extrapolate)
    max_progress_extrapolate = int(pm.max_progress + progress_extrapolate)
    progress_linspace_ex1 = np.linspace(min_progress_extrapolate, pm.min_progress, 20)
    progress_linspace_int = np.linspace(pm.min_progress, pm.max_progress, 60)
    progress_linspace_ex2 = np.linspace(pm.max_progress, max_progress_extrapolate, 20)

    # Calc min and max val in interval between 1% and 99% percentie
    min_val, max_val = pm.get_value_range([0.1, 0.9])
#     progress_linspace = np.linspace(min_progress_extrapolate, max_progress_extrapolate, 100)
#     min_val = float('inf')
#     max_val = float('-inf')
#     for quantile in [0.1, 0.9]:
#         curve = pm.get_quantile_curve(progress_linspace, quantile)
#         min_val = min(min_val, np.min(curve))
#         max_val = max(max_val, np.max(curve))

    #
    # Setup plot
    #
    biomarker_string = pt.get_biomarker_string(biomarker)
    figure_width = 6 if args.no_densities or args.only_densities else 12
    fig = plt.figure(figsize=(figure_width, 5))
    if args.only_densities:
        ax1 = None
        ax2 = plt.subplot(1, 1, 1)
        pt.setup_axes(plt, ax2, xgrid=False, ygrid=False)
    elif args.no_densities:
        ax1 = plt.subplot(1, 1, 1)
        ax2 = None
        pt.setup_axes(plt, ax1, xgrid=False, ygrid=False)
    else:
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
        pt.setup_axes(plt, ax1, xgrid=False, ygrid=False)
        pt.setup_axes(plt, ax2)

    if not args.only_densities:
        if args.no_model and not args.plot_synth_model:
            ax1.set_title('Aligned samples for {0}'.format(biomarker_string))
        else:
            ax1.set_title('Quantile curves for {0}'.format(biomarker_string))
        if args.phase == 'mciad':
            ax1.set_xlabel('Disease progress (days before/after conversion to AD)')
        else:
            ax1.set_xlabel('Disease progress (days before/after conversion to MCI)')
        ax1.set_ylabel(DataHandler.get_biomarker_unit(biomarker))
        if args.xlim is not None:
            ax1.set_xlim(args.xlim[0], args.xlim[1])
        else:
            ax1.set_xlim(min_progress_extrapolate, max_progress_extrapolate)
        if args.ylim is not None:
            ax1.set_ylim(args.ylim[0], args.ylim[1])

    #
    # Plot the percentile curves of the fitted model
    #
    if not args.no_model and not args.only_densities:
        ax1.axvline(pm.min_progress, color='0.15', linestyle=':')
        ax1.axvline(pm.max_progress, color='0.15', linestyle=':')

        quantiles = [0.1, 0.25, 0.5, 0.75, 0.9]
        grey_values = ['0.4', '0.2', '0', '0.2', '0.4']
        for grey_value, quantile in zip(grey_values, quantiles):
            curve_int = pm.get_quantile_curve(progress_linspace_int, quantile)
            ax1.plot(progress_linspace_int, curve_int, color=grey_value)

            if not args.no_extrapolation:
                curve_ex1 = pm.get_quantile_curve(progress_linspace_ex1, quantile)
                curve_ex2 = pm.get_quantile_curve(progress_linspace_ex2, quantile)
                ax1.plot(progress_linspace_ex1, curve_ex1, '--', color=grey_value)
                ax1.plot(progress_linspace_ex2, curve_ex2, '--', color=grey_value)

            if args.plot_quantile_label:
                label = '$q={0}\%$'.format(quantile * 100)
                ax1.text(progress_linspace_int[-1] + 10, curve_int[-1], label, fontsize=10)

        if args.plot_donohue:
            print 'Plotting Donohue'
            donohue_file = os.path.join(data_handler._conf.models_folder,
                                        'donohue', 'population_{0}.csv'.format(biomarker.replace(' ', '.')))
            if not os.path.isfile(donohue_file):
                print log.ERROR, 'Donohue model file not found: {0}'.format(donohue_file)
                return

            r = mlab.csv2rec(donohue_file)
            if args.method == 'joint':
                offset = 2200
            else:
                offset = 300
            progrs = r[r.dtype.names[0]] * 30.44 + offset
            vals = r[r.dtype.names[1]]
            curve_donohue = []
            progr_donohue = []
            for p in progress_linspace_int:
                if progrs[0] < p < progrs[-1]:
                    i = 1
                    while p > progrs[i]:
                        i += 1
                    # TODO linear interpolation
                    progr_donohue.append(progrs[i])
                    curve_donohue.append(vals[i])
            ax1.plot(progr_donohue, curve_donohue, '--', color='b', linewidth=2)

    #
    # Plot synthetic model curve
    #
    if plot_synth_model:
        progress_linspace_synth = np.linspace(-2500, 2500, 100)
        quantiles = [0.1, 0.25, 0.5, 0.75, 0.9]
        alphas = [0.4, 0.7, 1.0, 0.7, 0.4]
        for quantile, alpha in zip(quantiles, alphas):
            curve_synth = [SynthModel.get_distributed_value(biomarker, p, cdf=quantile) for p in progress_linspace_synth]
            ax1.plot(progress_linspace_synth, curve_synth, color='b', alpha=alpha)

    #
    # Plot predictor function
    #
    if args.plot_eta is not None and not args.only_densities:
        # Get second axis of plot 1
        ax1b = ax1.twinx()

        # Plot all progresses
        # ax1b.scatter(pm.all_progresses, pm.all_mus, facecolor='b', marker='o', edgecolor='none', alpha=0.2)
        ax1b.text(pm.progresses[-1], pm.sigmas[-1], '$\mu$', color='b', fontsize=11)

        # Plot binned progresses
        ax1b.scatter(pm.progresses, pm.sigmas, color='b', marker='x')

        # Plot interpolated model
        mus = [pm.get_eta(pm.sigmas, p) for p in progress_linspace_int]
        ax1b.plot(progress_linspace_int, mus, color='b')

        if not args.no_extrapolation:
            mus = [pm.get_eta(pm.sigmas, p) for p in progress_linspace_ex1]
            ax1b.plot(progress_linspace_ex1, mus, '--', color='b')
            mus = [pm.get_eta(pm.sigmas, p) for p in progress_linspace_ex2]
            ax1b.plot(progress_linspace_ex2, mus, '--', color='b')
        if args.xlim is not None:
            ax1b.set_xlim(args.xlim[0], args.xlim[1])
        else:
            ax1b.set_xlim(min_progress_extrapolate, max_progress_extrapolate)

    #
    # Plot errors
    #
    if args.plot_errors and not args.only_densities:
        eval_file = model_file.replace('.csv', '_eval_cover.csv')
        if not os.path.isfile(eval_file):
            print log.ERROR, 'Evaluation file not found: {0}'.format(eval_file)
        else:
            m = mlab.csv2rec(eval_file)
            progresses = m['progress']
            errors = m['error']

            # Get second axis of plot 1
            ax1b = ax1.twinx()
            # ax1b.set_ylim(0, max(150, 1.2 * np.max(errors)))
            ax1b.plot(progresses, errors, color='g', marker='x')
            ax1b.text(progresses[-1], errors[-1], 'Discr.', color='g', fontsize=11)
            ax1b.axhline(np.mean(errors), color='g', linestyle='--', alpha=0.5)

            median_curve = pm.get_quantile_curve(progresses, 0.5)
            min_value = np.min(median_curve)
            max_value = np.max(median_curve)
            rect = mpl.patches.Rectangle((progresses[0], min_value), progresses[-1] - progresses[0],
                                         max_value - min_value,
                                         fc=(0.0, 0.5, 0.0, 0.1), ec=(0.0, 0.5, 0.0, 0.8),
                                         linewidth=1)
            ax1.add_patch(rect)

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
            # diagn_points = [0.5 if p < 0 else 1.0 for p in progr_points]
            diagn_points = m['diagnosis']
            diagn_points[(0.25 <= diagn_points) & (diagn_points <= 0.75)] = 0.5

            print log.INFO, 'Plotting {0} sample points...'.format(len(progr_points))
            ax1.scatter(progr_points, value_points, s=15.0, c=diagn_points, edgecolor='none',
                        vmin=0.0, vmax=1.0, cmap=pt.progression_cmap, alpha=args.points_alpha)
            if args.phase == 'cnmci':
                rects = [mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_cn + (args.points_alpha,), linewidth=0),
                         mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_mci + (args.points_alpha,), linewidth=0)]
                labels = ['CN', 'MCI']
            elif args.phase == 'mciad':
                rects = [mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_mci + (args.points_alpha,), linewidth=0),
                         mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_ad + (args.points_alpha,), linewidth=0)]
                labels = ['MCI', 'AD']
            else:
                rects = [mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_cn + (args.points_alpha,), linewidth=0),
                         mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_mci + (args.points_alpha,), linewidth=0),
                         mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_ad + (args.points_alpha,), linewidth=0)]
                labels = ['CN', 'MCI', 'AD']
            legend = ax1.legend(rects, labels, fontsize=10, ncol=len(rects), loc='upper center', framealpha=0.9)
            legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    #
    # Plot PDFs
    #
    progr_samples = [-2000, -1000, 0, 1000, 2000, 3000, 4000] if args.phase == 'joint' else \
                    [-2000, -1500, -1000, -500, 0, 500, 1000, 1500, 2000]

    if args.phase == 'cnmci':
        vmin = -2000
        vmax = 6000
    elif args.phase == 'mciad':
        vmin = -6000
        vmax = 2000
    elif args.phase == 'joint':
        vmin = -2000
        vmax = 4000
    sample_cmap = cmx.ScalarMappable(
        norm=colors.Normalize(vmin=vmin, vmax=vmax),
        cmap=plt.get_cmap(pt.progression_cmap))

    if not args.no_sample_lines and not args.only_densities:
        for progr in progr_samples:
            if not args.no_extrapolation or pm.min_progress < progr < pm.max_progress:
                # sample_color = sample_cmap.to_rgba(progr_samples.index(progr))
                sample_color = sample_cmap.to_rgba(progr)
                linestyle = '--' if progr < pm.min_progress or progr > pm.max_progress else '-'
                ax1.axvline(progr, color=sample_color, linestyle=linestyle, alpha=0.3)

    if not args.no_densities:
        ax2.set_title('Probability density function for {0}'.format(biomarker_string))
        ax2.set_xlabel(DataHandler.get_biomarker_unit(biomarker))
        ax2.set_ylabel('Probability')
        if args.ylim is None:
            values = np.linspace(min_val, max_val, 250)
            ax2.set_xlim(min_val, max_val)
        else:
            values = np.linspace(args.ylim[0], args.ylim[1], 250)
            ax2.set_xlim(args.ylim[0], args.ylim[1])

        for progr in progr_samples:
            if not args.no_extrapolation or pm.min_progress < progr < pm.max_progress:
                # sample_color = sample_cmap.to_rgba(progr_samples.index(progr))
                sample_color = sample_cmap.to_rgba(progr)
                linestyle = '--' if progr < pm.min_progress or progr > pm.max_progress else '-'
                probs = pm.get_density_distribution(values, progr)
                ax2.plot(values, probs, label=str(progr), color=sample_color, linestyle=linestyle)

                if plot_synth_model:
                    probs = [SynthModel.get_probability(biomarker, progr, v) for v in values]
                    ax2.plot(values, probs, color='b', linestyle='--')

        legend = ax2.legend(fontsize=10, loc='best', framealpha=0.9)
        legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    #
    # Draw or save the plot
    #
    plt.tight_layout()
    if args.save_plots or args.plot_file is not None:
        if args.plot_file is not None:
            plot_filename = args.plot_file
        else:
            plot_filename = model_file.replace('.csv', '.pdf')
        plt.savefig(plot_filename, transparent=True)
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
