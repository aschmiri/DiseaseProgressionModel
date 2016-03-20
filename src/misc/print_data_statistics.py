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
    # print_gender_statistics(args, data_handler)
    print_terminal_decline_statistics(args, data_handler)


def print_terminal_decline_statistics(args, data_handler):
    # RIDS of subjects that withdrew from the study due to death.
    # These are all subjects with WDREASEON = 2 (ADNI1) or WDREASEON = 1 (ADNIGO/2) in TREATDIS.xlsx
    rids1 = {438, 103, 397, 1184, 884, 1338, 78, 1021, 1244, 825, 1277, 517, 821, 240, 1119, 177, 647, 67, 273, 786,
             559, 500, 607, 669, 293, 1211, 362, 963, 312, 1284, 57, 865, 155, 425, 326, 638, 1103}
    rids2 = {1203, 514, 4223, 4474, 15, 4237, 258, 289, 892, 830, 4609, 715, 408, 588, 4442, 4733, 376, 4770, 256, 294,
             108, 4892, 1271, 1394, 4282, 4897, 42, 1116, 4802, 1406, 1425, 947, 702, 4337, 4805, 649, 4910, 572, 388,
             4096, 1057, 922}
    ridsGO = {973, 1010, 1131, 1194, 2070, 128, 834, 845}

    # Subjects with death cause other than AD (Of the few where the death cause is actually indicated)
    rids_other_cause = {397, 78, 1021, 821, 647, 273, 963, 638,  # ADNI1
                        1203, 4892, 42, 4805,  # ADNI2
                        1131, 2070}  # ADNIGO

    rids_death_by_ad = rids1.union(rids2).union(ridsGO).difference(rids_other_cause)
    print '# of subjects that died due to AD or unknown cause: {0}'.format(len(rids_death_by_ad))

    measurements1 = data_handler.get_measurements_as_dict(min_visits=args.min_visits,
                                                          select_training_set=True,
                                                          no_regression=True)
    training_rids = set(measurements1.keys())
    print '# training subjects: {0}'.format(len(training_rids))
    print '  of which deceased: {0}'.format(len(training_rids.intersection(rids_death_by_ad)))

    measurements2 = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24'],
                                                          select_test_set=True,
                                                          select_complete=True,
                                                          no_regression=True)

    test_rids = set(measurements2.keys())
    print '# test subjects: {0}'.format(len(test_rids))
    print '  of which deceased: {0}'.format(len(test_rids.intersection(rids_death_by_ad)))

    assert len(training_rids.intersection(test_rids)) == 0

    # RIDS from RECADV.xlsx
    # rids2 = [884, 78, 517, 786, 1086, 293, 362, 638, 2119, 459, 834, 1407, 2106, 2379, 2233, 1066, 2087, 834, 1203,
    #          1271, 1080, 258, 4442, 1425, 256, 1331, 4802, 702]
    # rids2 = [1086, 2119, 459, 1407, 2106, 2379, 2233, 1066, 2087, 1080, 256, 1331]


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
