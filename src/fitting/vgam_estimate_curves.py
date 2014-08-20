#! /usr/bin/env python2.7
import os.path
import argparse
import csv
import subprocess
import joblib as jl
from common import log as log
from common import adni_tools as adni
from common import vgam as vgam


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser = vgam.add_common_arguments(parser)
    parser.add_argument('-n', '--nr_threads', type=int, default=4, help='number of threads')
    parser.add_argument('-d', '--degrees_of_freedom', type=int, default=2, help='degrees of freedom for the LMS method')
    parser.add_argument('--no_regression', action='store_true', default=False, help='do not perform age regression of biomarker values')
    args = parser.parse_args()

    # Get the data files and biomarkers
    data_files = vgam.get_data_files(args)
    biomarker_set = vgam.get_biomarker_set(args)

    # Estimate curves
    generate_csv_files(args, biomarker_set, data_files)
    estimate_model_all_biomarkers(args, biomarker_set)


def generate_csv_files(args, biomarkers, data_files):
    measurements = vgam.get_measurements_as_dict(data_files)
    data_folders = vgam.get_data_folders(args)

    if args.iteration > 0:
        input_folders = vgam.get_data_folders(args, previous_iteration=True)
        measurements = vgam.get_scaled_measurements(input_folders, measurements, biomarkers=['CDRSB'])

    for biomarker in biomarkers:
        print log.INFO, 'Generating output CSV for {0}...'.format(biomarker)
        data_folder = vgam.get_data_folder(data_folders, biomarker)
        csv_file = os.path.join(data_folder, biomarker.replace(' ', '_') + '.csv')
        writer = csv.writer(open(csv_file, 'wb'), delimiter=',')
        writer.writerow(['rid', 'progress', 'value'])

        subjects = set()
        num_samples = 0
        for rid, visits in measurements.items():
            for _, visit_data in visits.items():
                try:
                    progress = adni.safe_cast(visit_data['progress'], int)
                    biomarker_name = biomarker
                    value = adni.safe_cast(visit_data[biomarker_name], float)
                    if progress is not None and value is not None:
                        writer.writerow([rid, progress, value])
                        subjects.add(rid)
                        num_samples += 1
                except Exception:
                    pass

        print log.RESULT, 'Collected {0} samples from {1} subjects.'.format(num_samples, len(subjects))


def estimate_model_all_biomarkers(args, biomarkers):
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(estimate_model)(args, biomarker) for biomarker in biomarkers)


def estimate_model(args, biomarker):
    print log.INFO, 'Fitting curve to {0}...'.format(biomarker)
    r_file = os.path.join(adni.project_folder, 'src/fitting/vgam_estimate_curves.R')
    data_folders = vgam.get_data_folders(args)
    data_folder = vgam.get_data_folder(data_folders, biomarker)
    csv_file = os.path.join(data_folder, biomarker.replace(' ', '_') + '.csv')
    output_file = csv_file.replace('.csv', '_curves.csv')
    densities_file = csv_file.replace('.csv', '_densities.csv')
    image_file = csv_file.replace('.csv', '_curves.pdf')
    stdout_file = csv_file.replace('.csv', '_stdout.Rout')

    command = "R CMD BATCH \"--args input_file='%s' output_file='%s' degrees_of_freedom=%i save_plot=1 plot_file='%s' densities_file='%s'\" %s %s" \
              % (csv_file, output_file, int(args.degrees_of_freedom), image_file, densities_file, r_file, stdout_file)

    subprocess.call(command, shell=True)


if __name__ == '__main__':
    main()
