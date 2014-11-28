#! /usr/bin/env python2.7
import os.path
import argparse
import numpy as np
from sklearn import preprocessing
from sklearn import lda
from sklearn import svm
from sklearn import ensemble
from sklearn import cross_validation
from sklearn.metrics import make_scorer
from common import log as log
from common import evaluation_tools as et
from common.datahandler import DataHandler


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('visits', nargs='+', type=str, help='the viscodes to be sampled')
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-p', '--phase', default=None, choices=DataHandler.get_phase_choices(), help='the phase for which the model is to be trained')
    parser.add_argument('-c', '--classifier', default='svm', choices=['lda', 'svm', 'lsvm', 'rf'], help='the approach used to classify the subjects')
    parser.add_argument('--estimate_dprs', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--recompute_estimates', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--consistent_data', action='store_true', help='us only subjects with bl, m12 and m24 visits')
    parser.add_argument('--num_folds', type=int, default=10, help='number of folds for the n-fold cross validation')
    parser.add_argument('--latex_file', type=str, default=None, help='add output to a LaTeX file')
    args = parser.parse_args()

    # Get estimates
    rids, diagnoses, dpis, dprs, _, _ = et.get_progress_estimates(
        args.visits,
        method=args.method,
        biomarkers=args.biomarkers,
        phase=args.phase,
        estimate_dprs=args.estimate_dprs,
        recompute_estimates=args.recompute_estimates,
        consistent_data=args.consistent_data)

    # Select converters and non-converters sets
    rcds = get_rcds(args, rids, diagnoses, dpis, dprs)
    rfds = get_rfds(args, rids, diagnoses, dpis, dprs)

    # Analyse output
    analyse_decline(args, rids, dpis, dprs, rcds)
    analyse_decline(args, rids, dpis, dprs, rfds)


def get_rcds(args, rids, diagnoses, dpis, dprs):
    data_handler = DataHandler.get_data_handler()
    measurements = data_handler.get_measurements_as_dict(
        visits=['bl', 'm24'],
        biomarkers=['MMSE'],
        select_complete=True,
        no_regression=True)

    rcds = set()
    for rid in rids:
        if rid in measurements:
            mmse_bl = measurements[rid]['bl']['MMSE']
            mmse_m24 = measurements[rid]['m24']['MMSE']
            rcd = (mmse_bl - mmse_m24) >= 8
            if rcd:
                rcds.add(rid)

    print log.RESULT, 'Selected {0} subjects with rapid cognitive decline (RCD).'.format(len(rcds))
    return rcds


def get_rfds(args, rids, diagnoses, dpis, dprs):
    data_handler = DataHandler.get_data_handler()
    measurements = data_handler.get_measurements_as_dict(
        visits=['bl', 'm24'],
        biomarkers=['FAQ'],
        select_complete=True,
        no_regression=True)

    rfds = set()
    for rid in rids:
        if rid in measurements:
            faq_bl = measurements[rid]['bl']['FAQ']
            faq_m24 = measurements[rid]['m24']['FAQ']
            rcd = (faq_m24 - faq_bl) >= 10
            if rcd:
                rfds.add(rid)

    print log.RESULT, 'Selected {0} subjects with rapid functional decline (RFD).'.format(len(rfds))
    return rfds


def analyse_decline(args, rids, dpis, dprs, rds):
    print log.INFO, 'Analysing classification accuracies...'
    dpis = np.array(dpis)
    dprs = np.array(dprs)
    labels = np.array([1 if rid in rds else 0 for rid in rids])

    # Assemble features
    features = np.zeros((len(dpis), 2))
    features[:, 0] = dpis
    if args.estimate_dprs:
        features[:, 1] = dprs
    else:
        # Copy DPIs as second features as LDA needs two features
        features[:, 1] = dpis
    features = preprocessing.scale(features)

    acc, sens, spec = run_classification(args, features, labels)
    print log.RESULT, '{0}-fold cross validation, RD vs. non-RD ACC={1:.2f}, SENS={2:.2f}, SPEC={3:.2f}'.format(args.num_folds, acc, sens, spec)

    if args.latex_file is not None:
        data_handler = DataHandler.get_data_handler(method=args.method,
                                                    biomarkers=args.biomarkers,
                                                    phase=args.phase)
        filename = os.path.join(data_handler.get_eval_folder(), args.latex_file)
        print log.INFO, 'Writing classification results to {0}...'.format(filename)
        with open(filename, 'a') as latex_file:
            latex_file.write('{0} & {1} & {2:.2f} & {3:.2f} & {4:.2f}\\\\\n'.format(
                             args.method,
                             len(args.visits),
                             acc, sens, spec))


def sensitivity(ground_truth, predictions):
    TP = [1 if i == 1 and j == 1 else 0 for i, j in zip(ground_truth, predictions)]
    P = [1 if i == 1 else 0 for i in ground_truth]
    return float(np.sum(TP)) / np.sum(P)


def specificity(ground_truth, predictions):
    TN = [1 if i == 0 and j == 0 else 0 for i, j in zip(ground_truth, predictions)]
    N = [1 if i == 0 else 0 for i in ground_truth]
    return float(np.sum(TN)) / np.sum(N)


def run_classification(args, features, labels):
    class_weight = {0: 1.0, 1: float(len(labels) - sum(labels)) / float(sum(labels))}

    # Assemble data
    if args.classifier == 'lda':
        clf = lda.LDA()
    elif args.classifier == 'lsvm':
        clf = svm.LinearSVC()
    elif args.classifier == 'svm':
        clf = svm.SVC(kernel='rbf', class_weight=class_weight)
    elif args.classifier == 'rf':
        clf = ensemble.RandomForestClassifier()

    sensitivity_scorer = make_scorer(sensitivity, greater_is_better=True)
    specificity_scorer = make_scorer(specificity, greater_is_better=True)

    accuracies = cross_validation.cross_val_score(clf, features, labels, cv=args.num_folds, scoring='accuracy')
    sensitivities = cross_validation.cross_val_score(clf, features, labels, cv=args.num_folds, scoring=sensitivity_scorer)
    specificities = cross_validation.cross_val_score(clf, features, labels, cv=args.num_folds, scoring=specificity_scorer)

    return accuracies.mean(), sensitivities.mean(), specificities.mean()


if __name__ == '__main__':
    main()
