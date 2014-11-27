#! /usr/bin/env python2.7
import os.path
import argparse
import numpy as np
from sklearn import preprocessing
from sklearn import lda
from sklearn import svm
from sklearn import ensemble
from sklearn import cross_validation
from common import log as log
from common import evaluation_tools as et
from common.datahandler import DataHandler


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('visits', nargs='+', type=str, help='the viscodes to be sampled')
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-p', '--phase', default=None, choices=DataHandler.get_phase_choices(), help='the phase for which the model is to be trained')
    parser.add_argument('-c', '--classifier', default='lda', choices=['lda', 'svm', 'lsvm', 'rf'], help='the approach used to classify the subjects')
    parser.add_argument('--estimate_dprs', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--recompute_estimates', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--consistent_data', action='store_true', help='us only subjects with bl, m12 and m24 visits')
    parser.add_argument('--num_folds', type=int, default=10, help='number of folds for the n-fold cross validation')
    parser.add_argument('--latex_file', type=str, default=None, help='add output to a LaTeX file')
    args = parser.parse_args()

    # Get estimates
    _, diagnoses, dpis, dprs, _, _ = et.get_progress_estimates(
        args.visits,
        method=args.method,
        biomarkers=args.biomarkers,
        phase=args.phase,
        estimate_dprs=args.estimate_dprs,
        recompute_estimates=args.recompute_estimates,
        consistent_data=args.consistent_data)

    # Analyse estimates
    classify_diagnoses(args, dpis, dprs, diagnoses)


def classify_diagnoses(args, dpis, dprs, diagnoses):
    print log.INFO, 'Analysing classification accuracies...'

    # Calculate accuracy
    features = np.zeros((len(dpis), 2))
    features[:, 0] = np.array(dpis)
    if args.estimate_dprs:
        features[:, 1] = np.array(dprs)
    else:
        # Copy DPIs as second features as LDA needs two features
        features[:, 1] = np.array(dpis)
    features = preprocessing.scale(features)

    diagnoses = np.array(diagnoses)
    features_cn = features[diagnoses == 0.0]
    features_ad = features[diagnoses == 1.0]
    features_emci = features[diagnoses == 0.25]
    features_lmci = features[diagnoses == 0.75]
    features_mci = np.concatenate((features_emci, features_lmci))

    acc_cn_ad = run_classification(args, features_cn, features_ad)
    print log.RESULT, '{0}-fold cross validation, accuracy CN vs. AD {1}'.format(args.num_folds, acc_cn_ad)
    acc_cn_mci = run_classification(args, features_cn, features_mci)
    print log.RESULT, '{0}-fold cross validation, accuracy CN vs. MCI {1} '.format(args.num_folds, acc_cn_mci)
    acc_mci_ad = run_classification(args, features_mci, features_ad)
    print log.RESULT, '{0}-fold cross validation, accuracy MCI vs. AD {1}'.format(args.num_folds, acc_mci_ad)
    acc_emci_lmci = run_classification(args, features_emci, features_lmci)
    print log.RESULT, '{0}-fold cross validation, accuracy EMCI vs. LMCI {1}'.format(args.num_folds, acc_emci_lmci)
    acc_cn_mci_ad = run_classification(args, features_cn, features_mci, features_ad)
    print log.RESULT, '{0}-fold cross validation, accuracy CN vs. MCI vs. AD {1}'.format(args.num_folds, acc_cn_mci_ad)

    if args.latex_file is not None:
        data_handler = DataHandler.get_data_handler(method=args.method,
                                                    biomarkers=args.biomarkers,
                                                    phase=args.phase)
        filename = os.path.join(data_handler.get_eval_folder(), args.latex_file)
        print log.INFO, 'Writing classification results to {0}...'.format(filename)
        with open(filename, 'a') as latex_file:
            latex_file.write('{0} & {1} & {2:.2f} & {3:.2f} & {4:.2f} & {5:.2f} & {6:.2f}\\\\\n'.format(
                             args.method,
                             len(args.visits),
                             acc_cn_ad,
                             acc_cn_mci,
                             acc_mci_ad,
                             acc_emci_lmci,
                             acc_cn_mci_ad))


def run_classification(args, features1, features2, features3=None):
    # Assemble data
    if features3 is None:
        features = np.concatenate((features1, features2))
        labels = np.concatenate((np.zeros(len(features1)), np.ones(len(features2))))
    else:
        features = np.concatenate((features1, features2, features3))
        labels = np.concatenate((np.zeros(len(features1)), np.ones(len(features2)), 2 * np.ones(len(features3))))

    if args.classifier == 'lda':
        clf = lda.LDA()
    elif args.classifier == 'lsvm':
        clf = svm.LinearSVC()
    elif args.classifier == 'svm':
        clf = svm.SVC(kernel='rbf')
    elif args.classifier == 'rf':
        clf = ensemble.RandomForestClassifier()
    scores = cross_validation.cross_val_score(clf, features, labels, cv=args.num_folds)

    return scores.mean()  # , scores.std()


if __name__ == '__main__':
    main()
