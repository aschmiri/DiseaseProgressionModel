#!/usr/bin/env python
# print __doc__
import os.path
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cmx
from src.common import adni_tools as adni
from src.common import vgam as vgam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--biomarker_name', default=None, help='Name of the biomarker to be plotted')
    parser.add_argument('-p', '--no_points', action='store_true', default=False, help='Indicate that no points are to be plotted.')
    a = parser.parse_args()

    if a.biomarker_name is not None:
        biomarker_names = [a.biomarker_name]
    else:
        biomarker_names = adni.biomarker_names

    for biomarker in biomarker_names:
        print 'Generating plot for', biomarker

        points_file = os.path.join(adni.project_folder, 'data', biomarker.replace(' ', '_') + '.csv')
        curves_file = points_file.replace('.csv', '_curves.csv')
        if os.path.isfile(points_file) and os.path.isfile(curves_file):
            plot_model(biomarker, points_file, curves_file, not a.no_points)


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
        pfds = vgam.get_pfds_as_collection(biomarkers=[biomarker])
        y = pfds[biomarker]['values']

        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
        ax2.set_title('Probability density function for %s' % biomarker)
        ax2.set_xlabel('Volume' if biomarker in adni.volume_names else 'Score')
        ax2.set_ylabel('Probability')

        min_val = np.min(curves)
        max_val = np.max(curves)
        progr_samples = [-36, -18, 0, 18, 36]

        sample_cmap = cmx.ScalarMappable(
            norm=colors.Normalize(vmin=-len(progr_samples) + 1, vmax=(len(progr_samples) - 1)),
            cmap=plt.get_cmap(adni.adni_cmap))

        for progr in progr_samples:
            sample_color = sample_cmap.to_rgba(progr_samples.index(progr))
            ax1.axvline(progr, color=sample_color, linestyle='--', alpha=0.8)
            ax2.set_xlim(min_val, max_val)
            ax2.plot(y, pfds[biomarker][progr], label=str(progr), color=sample_color)

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
        ax1.scatter(progr_points, value_points, c=diagn_points, vmin=0.0, vmax=1.0, linewidths=0, cmap=adni.adni_cmap, alpha=0.25)

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
