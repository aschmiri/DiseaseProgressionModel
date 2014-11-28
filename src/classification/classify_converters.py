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
    _, _, dpis_conv, dprs_conv = select_converters(
        args, rids, diagnoses, dpis, dprs)
    _, _, dpis_nonconv, dprs_nonconv = select_nonconverters(
        args, rids, diagnoses, dpis, dprs)

    # Analyse data
    classify_converters(args, dpis_conv, dprs_conv, dpis_nonconv, dprs_nonconv)


def select_converters(args, rids, diagnoses, dpis, dprs):
    ''' Select data from subjects that convert within 2 years from MCI to AD. '''
    data_handler = DataHandler.get_data_handler(method=args.method)
    measurements = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24'],
                                                         no_regression=True,
                                                         select_training_set=True,
                                                         select_complete=True)

    # Select RIDSs of converters
    rids_select = set()
    for rid in measurements:
        if 0.25 <= measurements[rid]['bl']['DX.scan'] <= 0.75 and measurements[rid]['m24']['DX.scan'] == 1.0:
            rids_select.add(rid)

    selected_rids = []
    selected_diagnoses = []
    selected_dpis = []
    selected_dprs = []
    for i, rid in enumerate(rids):
        if rid in rids_select:
            selected_rids.append(rid)
            selected_diagnoses.append(diagnoses[i])
            selected_dpis.append(dpis[i])
            selected_dprs.append(dprs[i])

    print log.RESULT, 'Selected {0} converting subjects.'.format(len(selected_rids))
    return selected_rids, selected_diagnoses, selected_dpis, selected_dprs


def select_nonconverters(args, rids, diagnoses, dpis, dprs):
    ''' Select data from MCI subjects that do not convert. '''
    data_handler = DataHandler.get_data_handler(method=args.method)
    measurements = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24'],
                                                         no_regression=True,
                                                         select_test_set=True,
                                                         select_complete=True)
    # Select RIDSs of non-converters
    rids_select = set()
    for rid in measurements:
        if 0.25 <= measurements[rid]['bl']['DX.scan'] <= 0.75 and 0.25 <= measurements[rid]['m24']['DX.scan'] <= 0.75:
            rids_select.add(rid)

    selected_rids = []
    selected_diagnoses = []
    selected_dpis = []
    selected_dprs = []
    for i, rid in enumerate(rids):
        if rid in rids_select:
            selected_rids.append(rid)
            selected_diagnoses.append(diagnoses[i])
            selected_dpis.append(dpis[i])
            selected_dprs.append(dprs[i])

    print log.RESULT, 'Selected {0} non-converting subjects.'.format(len(selected_rids))
    return selected_rids, selected_diagnoses, selected_dpis, selected_dprs


def classify_converters(args, dpis_conv, dprs_conv, dpis_nonconv, dprs_nonconv):
    print log.INFO, 'Analysing classification accuracies...'
    dpis = np.concatenate((dpis_conv, dpis_nonconv))
    dprs = np.concatenate((dprs_conv, dprs_nonconv))
    labels = np.concatenate((np.ones(len(dpis_conv)), np.zeros(len(dpis_nonconv))))

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
    print log.RESULT, '{0}-fold cross validation, converters vs. non-converters ACC={1:.2f}, SENS={2:.2f}, SPEC={3:.2f}'.format(args.num_folds, acc, sens, spec)

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
