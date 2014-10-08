#! /usr/bin/env python2.7
import os.path
import pickle
from common import adni_tools as adni
from vgam.datahandler import DataHandler


def main():
    # Collect data for test
    data_handler = DataHandler.get_data_handler()
    biomarkers = data_handler.get_biomarker_set()

    mean_changes = {}
    for biomarker in biomarkers:
        measurements = data_handler.get_measurements_as_dict(visits=['bl', 'm12'],
                                                             biomarkers=[biomarker],
                                                             select_complete=True)

        mean_changes_biomarker = {0.0: 0.0, 0.25: 0.0, 0.5: 0.0, 0.75: 0.0, 1.0: 0.0}
        num_subjects = {0.0: 0, 0.25: 0, 0.5: 0, 0.75: 0, 1.0: 0}
        for rid in measurements:
            diagnosis = measurements[rid]['bl']['DX.scan']
            value_bl = measurements[rid]['bl'][biomarker]
            value_y1 = measurements[rid]['m12'][biomarker]
            scantime_bl = measurements[rid]['bl']['scantime']
            scantime_y1 = measurements[rid]['m12']['scantime']

            change = (value_y1 - value_bl) / (scantime_y1 - scantime_bl)

            mean_changes_biomarker[diagnosis] += change
            num_subjects[diagnosis] += 1

            if diagnosis in [0.25, 0.75]:
                mean_changes_biomarker[0.5] += change
                num_subjects[0.5] += 1

        for diagnosis in mean_changes_biomarker:
            mean_changes_biomarker[diagnosis] /= num_subjects[diagnosis]

        mean_changes.update({biomarker: mean_changes_biomarker})

    print mean_changes
    mean_changes_file = os.path.join(adni.eval_folder, 'mean_changes.p')
    pickle.dump(mean_changes, open(mean_changes_file, 'wb'))


if __name__ == '__main__':
    main()
