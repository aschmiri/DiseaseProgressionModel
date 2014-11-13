#! /usr/bin/env python2.7
import os.path
import argparse
import numpy as np
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
    parser.add_argument('--consistent_data', action='store_true', help='us only subjects with bl, m12 and m24 visits')
    parser.add_argument('--estimate_dprs', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--recompute_estimates', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--no_analysis', action='store_true', help='do not perform analysis')
    parser.add_argument('--no_plot', action='store_true', help='do not plot the results')
    parser.add_argument('--plot_distribution', action='store_true', help='plot the DP/DPR distribution')
    parser.add_argument('--plot_lines', action='store_true', help='plot graphs instead of matrix')
    parser.add_argument('--plot_steps', type=int, default=15, help='number of steps for the DPI scale')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output file')
    parser.add_argument('--plot_cmap_jet', action='store_true', help='use the colour map jet')
    parser.add_argument('--latex_file', type=str, default=None, help='add output to a LaTeX file')
    args = parser.parse_args()

    _, diagnoses, dpis, dprs, mean_min, mean_max = et.get_progress_estimates(args.visits,
                                                                             biomarkers=args.biomarkers,
                                                                             method=args.method,
                                                                             estimate_dprs=args.estimate_dprs,
                                                                             recompute_estimates=args.recompute_estimates,
                                                                             consistent_data=args.consistent_data)
    if not args.no_plot:
        plot_dpi_estimates(args, dpis, diagnoses, mean_min, mean_max)
    if args.plot_distribution and args.estimate_dprs:
        plot_dpi_dpr_distribution(args, dpis, dprs, diagnoses)
    if not args.no_analysis:
        analyse_dpi_estimates(args, dpis, dprs, diagnoses)


def plot_dpi_estimates(args, dpis, diagnoses, mean_min, mean_max):
    print log.INFO, 'Plotting estimates...'
    dpi_range = float(ModelFitter.TEST_DPI_MAX - ModelFitter.TEST_DPI_MIN)

    # Setup plot
    fig, ax = plt.subplots(figsize=(6, 2))
    biomarkers_str = args.method if args.biomarkers is None else ', '.join(args.biomarkers)
    ax.set_title('DP estimation using {0} at {1}'.format(biomarkers_str, ', '.join(args.visits)))
    ax.set_xlabel('DP')
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

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
    medians = []
    q25 = []
    q75 = []
    for diag in [0.0, 0.25, 0.75, 1.0]:
        row = diagnosis_indices[diag]
        num_subjects = np.sum(matrix[row])
        matrix[row] /= num_subjects

        indices = np.where(diagnoses == diag)
        median = np.median(dpis[indices])
        medians.append((median - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps)
        q25.append((median - np.percentile(dpis[indices], 25.0)) / dpi_range * args.plot_steps)
        q75.append((np.percentile(dpis[indices], 75.0) - median) / dpi_range * args.plot_steps)

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
    plt.axvline((mean_min - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps, color='k', linestyle=':', alpha=0.6)
    plt.axvline((mean_max - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps, color='k', linestyle=':', alpha=0.6)
    plt.axvline(-ModelFitter.TEST_DPI_MIN / dpi_range * args.plot_steps, color='k', linestyle='-', alpha=0.6)

    # Draw or save the plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def plot_dpi_dpr_distribution(args, dpis, dprs, diagnoses):
    print log.INFO, 'Plotting estimate distributions...'

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

    # Draw or save the plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def analyse_dpi_estimates(args, dpis, dprs, diagnoses):
    print log.INFO, 'Analysing classification accuracies...'

    # Calculate accuracy
    if args.estimate_dprs:
        features = np.zeros((len(dpis), 2))
        features[:, 0] = np.array(dpis)
        features[:, 1] = np.array(dprs)
    else:
        features = np.zeros((len(dpis), 1))
        features[:, 0] = np.array(dpis)
    diagnoses = np.array(diagnoses)

    features_cn = features[np.where(diagnoses == 0.0)]
    features_ad = features[np.where(diagnoses == 1.0)]
    features_emci = features[np.where(diagnoses == 0.25)]
    features_lmci = features[np.where(diagnoses == 0.75)]
    features_mci = np.concatenate((features_emci, features_lmci))

    # AD - CL
    acc_cn_ad, t_cn_ad = run_classification(features_cn, features_ad)
    print log.RESULT, 'Leave-one-out accuracy CN vs. AD {0} (threshold {1})'.format(acc_cn_ad, t_cn_ad)
    acc_cn_mci, t_cn_mci = run_classification(features_cn, features_mci)
    print log.RESULT, 'Leave-one-out accuracy CN vs. MCI {0} (threshold {1})'.format(acc_cn_mci, t_cn_mci)
    acc_mci_ad, t_mci_ad = run_classification(features_mci, features_ad)
    print log.RESULT, 'Leave-one-out accuracy MCI vs. AD {0} (threshold {1})'.format(acc_mci_ad, t_mci_ad)
    acc_emci_lmci, t_emci_lmci = run_classification(features_emci, features_lmci)
    print log.RESULT, 'Leave-one-out accuracy EMCI vs. LMCI {0} (threshold {1})'.format(acc_emci_lmci, t_emci_lmci)

    if args.latex_file is not None:
        filename = os.path.join(DataHandler.get_eval_folder(), args.latex_file)
        with open(filename, 'a') as latex_file:
            latex_file.write('{0} & {1} & {2:.2f} & {3:.1f} & {4:.2f} & {5:.1f} & {6:.2f} & {7:.1f} & {8:.2f} & {9:.1f}\\\\\n'.format(
                             args.method,
                             len(args.visits),
                             acc_cn_ad, t_cn_ad,
                             acc_cn_mci, t_cn_mci,
                             acc_mci_ad, t_mci_ad,
                             acc_emci_lmci, t_emci_lmci))


def run_classification(features1, features2):
    # Assemble data
    features = np.concatenate((features1, features2))
    labels = np.concatenate((np.zeros(len(features1)), np.ones(len(features2))))

#     from sklearn import cross_validation
#     from sklearn import svm
#     from sklearn import lda
#     from sklearn import tree
#
#     loo = cross_validation.LeaveOneOut(len(labels))
#
#     num_correct = 0
#     for train_index, test_index in loo:
#         X_train, X_test = features[train_index], features[test_index]
#         y_train, y_test = labels[train_index], labels[test_index]
# #         clf = svm.SVC()
# #         clf = lda.LDA()
#         clf = tree.DecisionTreeClassifier()
# #         clf = linear_model.LinearRegression()
#         clf.fit(X_train, y_train)
#
#         if clf.predict(X_test) == y_test:
#             num_correct += 1
#
#     mean_accuracy = float(num_correct) / len(labels)
#     mean_tresh = 0
#
#     return mean_accuracy, mean_tresh

    # Sort data
    indices = np.argsort(features[:, 0])
    features = features[indices]
    labels = labels[indices]

    num_correct = 0
    thresholds = []
    from sklearn import cross_validation
    loo = cross_validation.LeaveOneOut(len(features))

    # for test_i, test_label in enumerate(labels):
    for train_index, test_index in loo:
        X_train, X_test = features[train_index], features[test_index]
        y_train, y_test = labels[train_index], labels[test_index]

        # Train optimal threshold
        accuracies = []
        for thresh in range(len(y_train)):
            accuracies.append(thresh - np.sum(y_train[:thresh]) + np.sum(y_train[thresh:]))
        max_index = np.argmax(accuracies)

        # Get and store threshold
        threshold = 0.5 * (X_train[max_index - 1] + X_train[max_index])
        thresholds.append(threshold)

        # Test test sample
        if y_test == (0 if X_test <= threshold else 1):
            num_correct += 1

    mean_accuracy = float(num_correct) / len(labels)
    mean_tresh = np.mean(thresholds)

    return mean_accuracy, mean_tresh


if __name__ == '__main__':
    main()
