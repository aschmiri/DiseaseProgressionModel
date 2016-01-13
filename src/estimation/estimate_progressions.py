#! /usr/bin/env python2.7
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
from common import log as log
from common import evaluation_tools as et
from common import plotting_tools as pt
from common.datahandler import DataHandler
from common.modelfitter import ModelFitter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('visits', nargs='+', type=str, help='the viscodes to be sampled')
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-p', '--phase', default=None, choices=DataHandler.get_phase_choices(), help='the phase for which the model is to be trained')
    parser.add_argument('--estimate_dprs', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--recompute_estimates', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--consistent_data', action='store_true', help='us only subjects with bl, m12 and m24 visits')
    parser.add_argument('--no_plot', action='store_true', help='do not plot the results')
    parser.add_argument('--plot_lines', action='store_true', help='plot graphs instead of matrix')
    parser.add_argument('--plot_steps', type=int, default=15, help='number of steps for the DPI scale')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output file')
    parser.add_argument('--plot_cmap_jet', action='store_true', help='use the colour map jet')
    args = parser.parse_args()

    # Get estimates
    _, diagnoses, dpis, dprs, mean_min, mean_max = et.get_progress_estimates(
        args.visits,
        method=args.method,
        biomarkers=args.biomarkers,
        phase=args.phase,
        estimate_dprs=args.estimate_dprs,
        recompute_estimates=args.recompute_estimates,
        consistent_data=args.consistent_data)

    # Plot results
    if not args.no_plot:
        plot_dpi_estimates(args, dpis, diagnoses, mean_min, mean_max)
        if args.estimate_dprs:
            plot_dpi_dpr_distribution(args, dpis, dprs, diagnoses)


def plot_dpi_estimates(args, dpis, diagnoses, mean_min, mean_max):
    print log.INFO, 'Plotting estimates...'
    test_dpi_min, test_dpi_max, _ = ModelFitter.get_test_dpi_range(args.phase)
    dpi_range = float(test_dpi_max - test_dpi_min)
    dpi_factor = float(args.plot_steps) / dpi_range

    # Setup plot
    fig, ax = plt.subplots(figsize=(6, 2))
    biomarkers_str = args.method if args.biomarkers is None else ', '.join(args.biomarkers)
    ax.set_title('DP estimation using {0} at {1}'.format(biomarkers_str, ', '.join(args.visits)))
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    xticks = np.linspace(0, args.plot_steps, 7)
    ax.set_xticks(xticks)
    ax.set_xticklabels([int(float(tick) / dpi_factor + test_dpi_min) for tick in xticks])

    # Compute matrix
    diagnosis_indices = {0.0: 0, 0.25: 1, 0.5: 1, 0.75: 2, 1.0: 3}
    matrix = np.zeros((4, args.plot_steps + 1))
    for dpi, diag in zip(dpis, diagnoses):
        row = diagnosis_indices[diag]
        dpi_index = round((dpi - test_dpi_min) * dpi_factor)
        matrix[row, dpi_index] += 1.0

    # Draw annotations
    dpis = np.array(dpis)
    diagnoses = np.array(diagnoses)
    medians = []
    q25 = []
    q75 = []
    for diag in [0.0, 0.25, 0.75, 1.0]:
        row = diagnosis_indices[diag]
        matrix[row] /= np.sum(matrix[row])

        indices = np.where(diagnoses == diag)
        median = np.median(dpis[indices])
        medians.append((median - test_dpi_min) * dpi_factor)
        q25.append((median - np.percentile(dpis[indices], 25.0)) * dpi_factor)
        q75.append((np.percentile(dpis[indices], 75.0) - median) * dpi_factor)

    if args.plot_lines:
        ax.set_ylim(-0.01, 0.36)

        sample_cmap = cmx.ScalarMappable(
            norm=colors.Normalize(0.0, 1.0),
            cmap=plt.get_cmap(pt.progression_cmap))
        for diag in [0.0, 0.25, 0.75, 1.0]:
            row = diagnosis_indices[diag]
            plt.plot(matrix[row], color=sample_cmap.to_rgba(diag))
    else:
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels(['CN', 'EMCI', 'LMCI', 'AD'])

        cmap = plt.get_cmap('jet') if args.plot_cmap_jet else plt.get_cmap('Greys')
        barcol = 'w' if args.plot_cmap_jet else 'r'
        plt.errorbar(medians, [0, 1, 2, 3], xerr=[q25, q75], fmt=None, ecolor=barcol, elinewidth=2, capsize=4, capthick=2)
        plt.plot(medians, [0, 1, 2, 3], linestyle='', color=barcol, marker='|', markersize=15, markeredgewidth=2)
        plt.imshow(matrix, cmap=cmap, interpolation='nearest')
    plt.axvline((mean_min - test_dpi_min) * dpi_factor, color='k', linestyle=':', alpha=0.6)
    plt.axvline((mean_max - test_dpi_min) * dpi_factor, color='k', linestyle=':', alpha=0.6)
    plt.axvline((0.0 - test_dpi_min) * dpi_factor, color='k', linestyle='-', alpha=0.6)
    if args.phase == 'joint':
        data_handler = DataHandler.get_data_handler(method=args.method, biomarkers=args.biomarkers, phase=args.phase)
        plt.axvline((data_handler.get_model_offset() - test_dpi_min) * dpi_factor, color='k', linestyle='-', alpha=0.6)

    # Draw or save the plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def plot_dpi_dpr_distribution(args, dpis, dprs, diagnoses):
    print log.INFO, 'Plotting estimate distributions...'
    diagnoses = np.array(diagnoses)
    diagnoses[(0.25 <= diagnoses) & (diagnoses <= 0.75)] = 0.5

    # Setup plot
    fig, ax = plt.subplots()
    pt.setup_axes(plt, ax)

    biomarkers_str = args.method if args.biomarkers is None else ', '.join(args.biomarkers)
    ax.set_title('DP estimation using {0} at {1}'.format(biomarkers_str, ', '.join(args.visits)))
    ax.set_xlabel('DP')
    ax.set_ylabel('DPR')

    plt.scatter(dpis, dprs, c=diagnoses, edgecolor='none', s=25.0,
                vmin=0.0, vmax=1.0, cmap=pt.progression_cmap,
                alpha=0.5)

    # Plot legend
    rects = [mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_cn + (0.5,), linewidth=0),
             mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_mci + (0.5,), linewidth=0),
             mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_ad + (0.5,), linewidth=0)]
    labels = ['CN', 'MCI', 'AD']
    legend = ax.legend(rects, labels, fontsize=10, ncol=len(rects), loc='upper center', framealpha=0.9)
    legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    # Draw or save the plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
