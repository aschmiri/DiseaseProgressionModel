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
from common import vgam as vgam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('method', choices=['cog', 'reg', 'long', 'cons', 'graph', 'mbl'])
    parser.add_argument('-n', '--biomarker_name', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-p', '--no_points', action='store_true', default=False, help='indication that no points are to be plotted')
    args = parser.parse_args()

    if args.biomarker_name is not None:
        biomarker_names = [args.biomarker_name]
    else:
        biomarker_sets = {'cog': adni.cog_score_names,
                          'reg': adni.volume_names,
                          'long': adni.volume_names,
                          'cons': adni.volume_names,
                          'graph': adni.volume_names_essential,
                          'mbl': adni.manifold_coordinate_names}
        biomarker_names = biomarker_sets[args.method]

    for biomarker in biomarker_names:
        print log.INFO, 'Generating plot for {0}...'.format(biomarker)

        points_file = os.path.join(adni.project_folder, 'data', args.method, 'init', biomarker.replace(' ', '_') + '.csv')
        curves_file = points_file.replace('.csv', '_curves.csv')
        if os.path.isfile(points_file) and os.path.isfile(curves_file):
            plot_model(biomarker, points_file, curves_file, not args.no_points)


def plot_model(biomarker, points_file, curves_file, plot_points, save_file=False, plot_densities=True):
    plt.figure(figsize=(12, 5), dpi=100)

    #
    # Get data
    #
    r = mlab.csv2rec(curves_file)
    r_sorted = np.sort(r, order=r.dtype.names[2])
    progrs = r_sorted[r.dtype.names[2]]
    curves = [r_sorted['x1'], r_sorted['x5'], r_sorted['x25'], r_sorted['x50'],
              r_sorted['x75'], r_sorted['x95'], r_sorted['x99']]

    #
    # Plot PDFs
    #
    if plot_densities is not None:
        pdfs = vgam.get_pfds_as_collection(folder=os.path.dirname(curves_file),
                                           biomarkers=[biomarker])

        values = pdfs[biomarker].pop('values')

        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
        ax2.set_title('Probability density function for %s' % biomarker)
        ax2.set_xlabel('Volume' if biomarker in adni.volume_names else 'Score')
        ax2.set_ylabel('Probability')

        min_val = np.min(curves)
        max_val = np.max(curves)
        progr_samples = [-1100, -550, 0, 500, 1100]

        sample_cmap = cmx.ScalarMappable(
            norm=colors.Normalize(vmin=-len(progr_samples) + 1, vmax=(len(progr_samples) - 1)),
            cmap=plt.get_cmap(aplt.progression_cmap))

        for i, progr in enumerate(progr_samples):
            pdf = vgam.interpolate_pdf(pdfs[biomarker], progr)
            sample_color = sample_cmap.to_rgba(i)
            ax1.axvline(progr, color=sample_color, linestyle='--', alpha=0.8)
            ax2.set_xlim(min_val, max_val)
            ax2.plot(values, pdf, label=str(progr), color=sample_color)

        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(handles2, labels2, fontsize=10)
    else:
        ax1 = plt.subplot(1, 1, 1)

    #
    # Plot points
    #
    if plot_points:
        m = mlab.csv2rec(points_file)
        progr_points = m['progress']
        value_points = m['value']
        diagn_points = [0.5 if p < 0 else 1.0 for p in progr_points]

        print log.INFO, 'Plotting {0} sample points...'.format(len(progr_points))
        ax1.scatter(progr_points, value_points, c=diagn_points, vmin=0.0, vmax=1.0, linewidths=0, cmap=aplt.progression_cmap, alpha=0.25)

    #
    # Plot the percentile curves
    #
    labels = ['1%', '5%', '25%', '50%', '75%', '95%', '99%']
    greyvals = ['0.45', '0.3', '0.15', '0', '0.15', '0.3', '0.45']
    # styles   = ['g-', 'g-', 'g-', 'b-', 'g-', 'g-', 'g-']

    for (curve, greyval, label) in zip(curves, greyvals, labels):
        ax1.plot(progrs, curve, color=greyval)
        ax1.text(progrs[-1] + 5, curve[-1] - 0.5, label, fontsize=10)

    ax1.set_title('Percentile curves for %s' % biomarker)
    ax1.set_xlabel('Disease progression relative to point of conversion')
    ax1.set_ylabel('Volume' if biomarker in adni.volume_names else 'Score')
    progress_offset = (vgam.MAX_PROGRESS - vgam.MIN_PROGRESS) * 0.1
    ax1.set_xlim(vgam.MIN_PROGRESS - progress_offset, vgam.MAX_PROGRESS + 2 * progress_offset)
    ax1.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.8, 0.8, 0.0), linewidth=0),
                mpl.patches.Rectangle((0, 0), 1, 1, fc=(1.0, 0.0, 0.0), linewidth=0)],
               ['MCI', 'AD'], fontsize=10)

    plt.tight_layout()

    #
    # Draw or save the plot
    #
    if save_file:
        plot_filename = curves_file.replace('.csv', '_plot.png')
        plt.savefig(plot_filename, dpi=100)
    else:
        plt.show()


if __name__ == '__main__':
    main()


# def plot_model_analytic(biomarker, points_file, curves_file, plot_points,
#                         save_file=False, plot_densities=True, plot_mu=True):
#     plt.figure(figsize=(12, 5), dpi=100)
#
#     #
#     # Get data
#     #
#     r = mlab.csv2rec(curves_file)
#     r_sorted = np.sort(r, order=r.dtype.names[2])
#     progrs = r_sorted[r.dtype.names[2]]
#     curves = [r_sorted['x1'], r_sorted['x5'], r_sorted['x25'], r_sorted['x50'],
#               r_sorted['x75'], r_sorted['x95'], r_sorted['x99']]
#
#     parameters = {}
#     for i in range(len(progrs)):
#         progr = progrs[i]
#         if progr not in parameters:
#             parameters.update({progr: {}})
#             parameters[progr].update({'sigma': np.exp(r_sorted['logsigma'][i])})
#             parameters[progr].update({'lambda': r_sorted['lambda'][i]})
#             parameters[progr].update({'mu': r_sorted['mu'][i]})
#             parameters[progr].update({'yoffset': r_sorted['fitmiscyoffset'][i]})
#
#     #
#     # Plot PDFs
#     #
#     if plot_densities is not None:
#         ax1 = plt.subplot(1, 2, 1)
#         ax2 = plt.subplot(1, 2, 2)
#         ax2.set_title('Probability density function for %s' % biomarker)
#         ax2.set_xlabel('Volume' if biomarker in adni.volume_names else 'Score')
#         ax2.set_ylabel('Probability')
#
#         min_val = np.min(curves)
#         max_val = np.max(curves)
#         values = np.arange(min_val, max_val, 100)
#         progr_samples = [-36, -18, 3, 18, 33]
#
#         sample_cmap = cmx.ScalarMappable(
#             norm=colors.Normalize(vmin=-len(progr_samples) + 1, vmax=(len(progr_samples) - 1)),
#             cmap=plt.get_cmap(adni.adni_cmap))
#
#         for progr in progr_samples:
#             sample_color = sample_cmap.to_rgba(progr_samples.index(progr))
#             ax1.axvline(progr, color=sample_color, linestyle='--', alpha=0.8)
#             ax2.set_xlim(min_val, max_val)
#             probs = [vgam.yeojohnson_density(parameters[progr]['yoffset'] + v,
#                                              parameters[progr]['lambda'],
#                                              parameters[progr]['mu'],
#                                              parameters[progr]['sigma']) for v in values]
#             ax2.plot(values, probs, label=str(progr), color=sample_color)
#
#         handles2, labels2 = ax2.get_legend_handles_labels()
#         ax2.legend(handles2, labels2, fontsize=10)
#     else:
#         ax1 = plt.subplot(1, 1, 1)
#
#     #
#     # Plot points
#     #
#     if plot_points:
#         m = mlab.csv2rec(points_file)
#         progr_points = m['progress']
#         value_points = m['value']
#         diagn_points = [0.5 if p < 0 else 1.0 for p in progr_points]
#         ax1.scatter(progr_points, value_points, c=diagn_points, vmin=0.0, vmax=1.0, linewidths=0, cmap=adni.adni_cmap, alpha=0.25)
#
#     #
#     # Plot mu
#     #
#     if plot_mu:
#         m = [parameters[progr]['mu'] for progr in progrs]
#         ax1.plot(progrs, m, color='blue')
#
#     #
#     # Plot the percentile curves
#     #
#     labels = ['1%', '5%', '25%', '50%', '75%', '95%', '99%']
#     greyvals = ['0.45', '0.3', '0.15', '0', '0.15', '0.3', '0.45']
#     # styles   = ['g-', 'g-', 'g-', 'b-', 'g-', 'g-', 'g-']
#
#     for (curve, greyval, label) in zip(curves, greyvals, labels):
#         ax1.plot(progrs, curve, color=greyval)
#         ax1.text(progrs[-1] + 5, curve[-1] - 0.5, label, fontsize=10)
#
#     ax1.set_title('Percentile curves for %s' % biomarker)
#     ax1.set_xlabel('Disease progression relative to point of conversion')
#     ax1.set_ylabel('Volume' if biomarker in adni.volume_names else 'Score')
#     progress_offset = (vgam.MAX_PROGRESS - vgam.MIN_PROGRESS) * 0.1
#     ax1.set_xlim(vgam.MIN_PROGRESS - progress_offset, vgam.MAX_PROGRESS + 2 * progress_offset)
#     ax1.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.8, 0.8, 0.0), linewidth=0),
#                 mpl.patches.Rectangle((0, 0), 1, 1, fc=(1.0, 0.0, 0.0), linewidth=0)],
#                ['MCI', 'AD'], fontsize=10)
#
#     plt.tight_layout()
#
#     #
#     # Draw or save the plot
#     #
#     if save_file:
#         plot_filename = curves_file.replace('.csv', '_plot.png')
#         plt.savefig(plot_filename, dpi=100)
#     else:
#         plt.show()
