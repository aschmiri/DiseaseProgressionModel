#! /usr/bin/env python2.7
import os.path
import argparse
import csv
import numpy as np
from common import log as log
from common.datahandler import DataHandler


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-p', '--phase', default=None, choices=DataHandler.get_phase_choices(), help='the phase for which the model is to be trained')
    parser.add_argument('-n', '--nr_threads', type=int, default=1, help='number of threads')
    parser.add_argument('--min_visits', type=int, default=0, help='the minimal number of visits')
    parser.add_argument('--no_regression', action='store_true', default=False, help='do not perform age regression of biomarker values')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    args = parser.parse_args()

    # Get the data files and biomarkers
    data_handler = DataHandler.get_data_handler(method=args.method,
                                                biomarkers=args.biomarkers,
                                                phase=args.phase)

    # Estimate curves
    # generate_csv_file(args, data_handler)
    print_gender_statistics(args, data_handler)


def print_visits_statistics(args, data_handler):
    measurements = data_handler.get_measurements_as_dict(min_visits=args.min_visits,
                                                         select_training_set=True,
                                                         no_regression=True)

    num_visits = []
    for rid, visits in measurements.items():
        num_visits.append(len(visits.items()))

    print log.RESULT, 'Min number of visits: {0}'.format(np.min(num_visits))
    print log.RESULT, 'Max number of visits: {0}'.format(np.max(num_visits))
    print log.RESULT, 'Mean number of visits: {0}'.format(np.median(num_visits))
    print log.RESULT, 'Median number of visits: {0}'.format(np.mean(num_visits))


def print_training_samples_statistics(args, data_handler):
    biomarkers = data_handler.get_biomarker_names()
    measurements = data_handler.get_measurements_as_dict(min_visits=args.min_visits,
                                                         select_training_set=True,
                                                         no_regression=True)
    for biomarker in biomarkers:
        subjects = set()
        num_samples = 0
        for rid, visits in measurements.items():
            for _, visit_data in visits.items():
                try:
                    progress = DataHandler.safe_cast(visit_data['progress'], int)
                    value = DataHandler.safe_cast(visit_data[biomarker], float)
                    if progress is not None and value is not None:
                        subjects.add(rid)
                        num_samples += 1
                except KeyError:
                    pass

        print log.RESULT, 'Biomarker {0}: collected {1} samples from {2} subjects.'.format(biomarker, num_samples, len(subjects))


def print_gender_statistics(args, data_handler):
    measurements = data_handler.get_measurements_as_dict(min_visits=args.min_visits,
                                                         select_training_set=True,
                                                         no_regression=True)

    # Setup DB
    import rpy2.robjects as robjects
    robjects.r['load'](os.path.join(data_handler.get_data_folder(), 'adnimerge/adnimerge.rdata'))
    adnimerge = robjects.r['adnimerge']

    # Print each row
    dx_dict = {0.0: 'CN', 0.25: 'EMCI', 0.5: 'MCI', 0.75: 'LMCI', 1.0: 'AD'}
    males = {'CN': 0, 'EMCI': 0, 'LMCI': 0, 'AD': 0}
    females = {'CN': 0, 'EMCI': 0, 'LMCI': 0, 'AD': 0}
    rids = adnimerge[adnimerge.colnames.index('RID')]
    for rid in measurements:
        gender = None
        for row in range(adnimerge.nrow):
            if rids[row] == rid:
                gender = adnimerge[9][row]
                dx = dx_dict[measurements[rid]['bl']['DX.scan']]
                break
        if gender is None:
            print log.WARNING, 'Gender could not be determined for RID {0}'.format(rid)

        if gender == 1:
            males[dx] += 1
        else:
            females[dx] += 1

    print log.RESULT, 'Training data:'
    print log.RESULT, '  Males:   {0} ({1})'.format(np.sum(males.values()), males)
    print log.RESULT, '  Females: {0} ({1})'.format(np.sum(females.values()), females)

    measurements = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24'],
                                                         select_test_set=True,
                                                         select_complete=True,
                                                         no_regression=True)

    # Print each row
    dx_dict = {0.0: 'CN', 0.25: 'EMCI', 0.5: 'MCI', 0.75: 'LMCI', 1.0: 'AD'}
    males = {'CN': 0, 'EMCI': 0, 'LMCI': 0, 'AD': 0}
    females = {'CN': 0, 'EMCI': 0, 'LMCI': 0, 'AD': 0}
    rids = adnimerge[adnimerge.colnames.index('RID')]
    for rid in measurements:
        if 'bl' not in measurements[rid]:
            print log.WARNING, 'BL not available for RID {0}'.format(rid)
            print measurements[rid]
            continue

        gender = None
        for row in range(adnimerge.nrow):
            if rids[row] == rid:
                gender = adnimerge[9][row]
                dx = dx_dict[measurements[rid]['bl']['DX.scan']]
                break
        if gender is None:
            print log.WARNING, 'Gender could not be determined for RID {0}'.format(rid)

        if gender == 1:
            males[dx] += 1
        else:
            females[dx] += 1

    print log.RESULT, 'Test data:'
    print log.RESULT, '  Males:   {0} ({1})'.format(np.sum(males.values()), males)
    print log.RESULT, '  Females: {0} ({1})'.format(np.sum(females.values()), females)


def write_biomarker_data(args, data_handler):
    dx_dict = {0.0: 'CN', 0.25: 'EMCI', 0.5: 'MCI', 0.75: 'LMCI', 1.0: 'AD'}
    biomarkers = data_handler.get_biomarker_names()

    # Print test data
    measurements = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24'],
                                                         select_test_set=True,
                                                         no_regression=True,
                                                         select_complete=True)
    print 'Selected {0} test subjects.'.format(len(measurements))

    with open('../biomarkers/volumes2.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        row = ['RID', 'VISCODE', 'DX.scan', 'AGE.scan', 'ScanDate']
        row += biomarkers
        spamwriter.writerow(row)

        for rid in measurements:
            for viscode in ['bl', 'm12', 'm24']:
                row = [rid,
                       viscode,
                       dx_dict[measurements[rid][viscode]['DX.scan']],
                       measurements[rid][viscode]['AGE.scan'],
                       measurements[rid][viscode]['scandate']]

                for biomarker in biomarkers:
                    if biomarker in measurements[rid][viscode]:
                        row.append(measurements[rid][viscode][biomarker])
                    else:
                        row.append('')

                spamwriter.writerow(row)

    # Print training data
    measurements = data_handler.get_measurements_as_dict(biomarkers=biomarkers,
                                                         select_training_set=True,
                                                         no_regression=True)
    print 'Selected {0} training subjects.'.format(len(measurements))

    with open('../biomarkers/volumes3.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        row = ['RID', 'VISCODE', 'DX.scan', 'AGE.scan', 'ScanDate']
        row += biomarkers
        spamwriter.writerow(row)

        for rid in measurements:
            for viscode in measurements[rid]:
                row = [rid,
                       viscode,
                       dx_dict[measurements[rid][viscode]['DX.scan']],
                       measurements[rid][viscode]['AGE.scan'],
                       measurements[rid][viscode]['scandate']]

                for biomarker in biomarkers:
                    if biomarker in measurements[rid][viscode]:
                        row.append(measurements[rid][viscode][biomarker])
                    else:
                        row.append('')

                spamwriter.writerow(row)


if __name__ == '__main__':
    main()
