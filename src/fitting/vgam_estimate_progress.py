#! /usr/bin/env python2.7
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt
from common import log as log
from common import adni_tools as adni
from vgam.datahandler import DataHandler
from vgam.modelfitter import ModelFitter
import fitting.vgam_evaluation as ve


def main():
    parser = argparse.ArgumentParser()
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('visits', nargs='+', type=str, help='the viscodes to be sampled')
    parser.add_argument('--consistent_data', action='store_true', help='us only subjects with bl, m12 and m24 visits')
    parser.add_argument('--estimate_dprs', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--recompute_estimates', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--no_plot', action='store_true', help='do not plot the results')
    parser.add_argument('--plot_steps', type=int, default=15, help='number of steps for the DPI scale')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output file')
    parser.add_argument('--latex_file', type=str, default=None, help='add output to a LaTeX file')
    args = parser.parse_args()

    _, diagnoses, dpis, _, mean_min, mean_max = ve.get_progress_estimates(args, args.visits,
                                                                          estimate_dprs=args.estimate_dprs,
                                                                          recompute_estimates=args.recompute_estimates,
                                                                          consistent_data=args.consistent_data)
    if not args.no_plot:
        plot_dpi_estimates(args, dpis, diagnoses, mean_min, mean_max)
    analyse_dpi_estimates(args, dpis, diagnoses)


def plot_dpi_estimates(args, dpis, diagnoses, mean_min, mean_max):
    dpi_range = float(ModelFitter.TEST_DPI_MAX - ModelFitter.TEST_DPI_MIN)

    # Setup plot
    fig, ax = plt.subplots(figsize=(6, 2))
    biomarkers_str = args.method if args.biomarkers is None else ', '.join(args.biomarkers)
    ax.set_title('DPI estimation using {0} at {1}'.format(biomarkers_str, ', '.join(args.visits)))
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
    for dpi, diag in zip(dpis, diagnoses):
        row = diagnosis_indices[diag]
        dpi_index = round((dpi - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps)
        matrix[row, dpi_index] += 1.0

    # Draw annotations
    dpis = np.array(dpis)
    diagnoses = np.array(diagnoses)
    means = []
    stds = []
    for diag in [0.0, 0.25, 0.75, 1.0]:
        row = diagnosis_indices[diag]
        num_subjects = np.sum(matrix[row])
        matrix[row] /= num_subjects

        indices = np.where(diagnoses == diag)
        means.append((np.mean(dpis[indices]) - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps)
        stds.append(np.std(dpis[indices]) / dpi_range * args.plot_steps)
        # ax.text(args.plot_steps + 1, i, '$n={0}$'.format(int(num_subjects)))

    plt.errorbar(means, [0, 1, 2, 3], xerr=stds, fmt=None, ecolor='r', elinewidth=2, capsize=4, capthick=2)
    plt.plot(means, [0, 1, 2, 3], linestyle='', color='r', marker='|', markersize=15, markeredgewidth=2)
    plt.axvline((mean_min - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps, color='k', linestyle=':', alpha=0.6)
    plt.axvline((mean_max - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps, color='k', linestyle=':', alpha=0.6)
    plt.axvline(-ModelFitter.TEST_DPI_MIN / dpi_range * args.plot_steps, color='k', linestyle='-', alpha=0.6)
    plt.imshow(matrix, cmap=plt.get_cmap('Greys'), interpolation='nearest')

    # Draw or save the plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def analyse_dpi_estimates(args, dpis, diagnoses):
    # Calculate accuracy
    dpis = np.array(dpis)
    diagnoses = np.array(diagnoses)

    dpis_cn = dpis[np.where(diagnoses == 0.0)]
    dpis_ad = dpis[np.where(diagnoses == 1.0)]
    dpis_emci = dpis[np.where(diagnoses == 0.25)]
    dpis_lmci = dpis[np.where(diagnoses == 0.75)]
    dpis_mci = np.concatenate((dpis_emci, dpis_lmci))

    # AD - CL
    acc_cn_ad, t_cn_ad = run_classification(dpis_cn, dpis_ad)
    print log.RESULT, 'Leave-one-out accuracy CN vs. AD {0} (threshold {1})'.format(acc_cn_ad, t_cn_ad)
    acc_cn_mci, t_cn_mci = run_classification(dpis_cn, dpis_mci)
    print log.RESULT, 'Leave-one-out accuracy CN vs. MCI {0} (threshold {1})'.format(acc_cn_mci, t_cn_mci)
    acc_mci_ad, t_mci_ad = run_classification(dpis_mci, dpis_ad)
    print log.RESULT, 'Leave-one-out accuracy MCI vs. AD {0} (threshold {1})'.format(acc_mci_ad, t_mci_ad)
    acc_emci_lmci, t_emci_lmci = run_classification(dpis_emci, dpis_lmci)
    print log.RESULT, 'Leave-one-out accuracy EMCI vs. LMCI {0} (threshold {1})'.format(acc_emci_lmci, t_emci_lmci)

    if args.latex_file is not None:
        filename = os.path.join(adni.eval_folder, args.latex_file)
        with open(filename, 'a') as latex_file:
            latex_file.write('{0} & {1} & {2:.2f} & {3:.1f} & {4:.2f} & {5:.1f} & {6:.2f} & {7:.1f} & {8:.2f} & {9:.1f}\\\\\n'.format(
                             args.method,
                             len(args.visits),
                             acc_cn_ad, t_cn_ad,
                             acc_cn_mci, t_cn_mci,
                             acc_mci_ad, t_mci_ad,
                             acc_emci_lmci, t_emci_lmci)
)

def run_classification(dpis1, dpis2):
    # Assemble data
    dpis = np.concatenate((dpis1, dpis2))
    labels = np.concatenate((np.zeros(len(dpis1)), np.ones(len(dpis2))))

    # Sort data
    indices = np.argsort(dpis)
    dpis = dpis[indices]
    labels = labels[indices]

    num_correct = 0
    thresholds = []
    for test_i, test_label in enumerate(labels):
        # Get training data
        train_dpis = np.delete(dpis, test_i)
        train_labels = np.delete(labels, test_i)

        # Train optimal threshold
        accuracies = []
        for thresh in range(len(train_labels)):
            accuracies.append(thresh - np.sum(train_labels[:thresh]) + np.sum(train_labels[thresh:]))
        max_index = np.argmax(accuracies)

        # Get and store threshold
        threshold = 0.5 * (train_dpis[max_index - 1] + train_dpis[max_index])
        thresholds.append(threshold)

        # Test test sample
        if test_label == (0 if dpis[test_i] < threshold else 1):
            num_correct += 1

    mean_accuracy = float(num_correct) / len(labels)
    mean_tresh = np.mean(thresholds)

    return mean_accuracy, mean_tresh


if __name__ == '__main__':
    main()
