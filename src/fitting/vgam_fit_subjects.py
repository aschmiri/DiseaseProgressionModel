#! /usr/bin/env python2.7
import argparse
import numpy as np
import matplotlib.pyplot as plt
from common import log as log
from vgam.datahandler import DataHandler
from vgam.modelfitter import ModelFitter
import fitting.vgam_evaluation as ve


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('visits', nargs='+', type=str, help='the viscodes to be sampled')
    parser.add_argument('--recompute_dpis', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--plot_steps', type=int, default=15, help='number of steps for the DPI scale')
    parser.add_argument('--output_file', type=str, default=None, help='filename of the output file')
    args = parser.parse_args()

    _, diagnoses, dpis, mean_min, mean_max = ve.get_dpi_estimates(args)
    plot_dpi_estimates(args, dpis, diagnoses, mean_min, mean_max)
    analyse_dpi_estimates(args, dpis, diagnoses)


def plot_dpi_estimates(args, dpis, diagnoses, mean_min, mean_max):
    dpi_range = float(ModelFitter.TEST_DPI_MAX - ModelFitter.TEST_DPI_MIN)

    # Setup plot
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.set_title('DPI estimation using {0} at {1}'.format(args.method, ', '.join(args.visits)))
    ax.set_xlabel('DPI')
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(['CN', 'EMCI', 'LMCI', 'AD'])

    xticks = np.linspace(0, args.plot_steps, 7)
    ax.set_xticks(xticks)
    ax.set_xticklabels([int(float(tick) / args.plot_steps * dpi_range + ModelFitter.TEST_DPI_MIN) for tick in xticks])

    # Compute matrix
    diagnosis_indices = {0.0: 0, 0.25: 1, 0.5: 1, 0.75: 2, 1.0: 3}
    matrix = np.zeros((4, args.plot_steps + 1))
    mean = np.zeros(4)
    for dpi, diag in zip(dpis, diagnoses):
        diagnosis_index = diagnosis_indices[diag]
        dpi_index = round((dpi - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps)
        matrix[diagnosis_index, dpi_index] += 1.0
        mean[diagnosis_index] += dpi

    # Draw annotations
    dpis = np.array(dpis)
    diagnoses = np.array(diagnoses)
    means = []
    stds = []
    for diag in [0.0, 0.25, 0.75, 1.0]:
        i = diagnosis_indices[diag]
        num_subjects = np.sum(matrix[i])
        matrix[i] /= num_subjects

        indices = np.where(diagnoses == diag)
        mean = (np.mean(dpis[indices]) - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps
        std = np.std(dpis[indices]) / dpi_range * args.plot_steps
        means.append(mean)
        stds.append(std)
        # ax.text(args.plot_steps + 1, i, '$n={0}$'.format(int(num_subjects)))

    plt.errorbar(means, [0, 1, 2, 3], xerr=stds, fmt=None, ecolor='r', elinewidth=2, capsize=4, capthick=2)
    plt.plot(means, [0, 1, 2, 3], linestyle='', color='r', marker='|', markersize=15, markeredgewidth=2)
    plt.axvline((mean_min - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps, color='k', linestyle=':', alpha=0.6)
    plt.axvline((mean_max - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps, color='k', linestyle=':', alpha=0.6)
    plt.axvline(-ModelFitter.TEST_DPI_MIN / dpi_range * args.plot_steps, color='k', linestyle='-', alpha=0.6)
    plt.imshow(matrix, cmap=plt.get_cmap('Greys'), interpolation='nearest')

    # Draw or save the plot
    plt.tight_layout()
    if args.output_file is not None:
        plt.savefig(args.output_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def analyse_dpi_estimates(args, dpis, diagnoses):
    # Calculate accuracy
    misclassified = 0
    for dpi, diagnosis in zip(dpis, diagnoses):
        if dpi > 0 and diagnosis == 0.5 or dpi <= 0 and diagnosis == 1.0:
            misclassified += 1
    print log.RESULT, 'Accuracy: {0}'.format(1.0 - float(misclassified) / float(len(dpis)))


if __name__ == '__main__':
    main()
