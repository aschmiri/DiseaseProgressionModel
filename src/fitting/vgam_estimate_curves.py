#!/usr/bin/env python

import argparse
import os.path
import csv
import subprocess
import joblib as jl
from src.common import adni_tools as adni
from src.common import vgam as vgam


def main():
    # parse command line options
    parser = argparse.ArgumentParser(description='Generate reference curves and reports for BAMBI.')
    parser.add_argument('-n', '--nr_threads', dest='nr_threads', type=int, default=4, help='number of threads')
    parser.add_argument('-d', '--degrees-of-freedom', dest='degrees_of_freedom', type=int, default=2, help='degrees of freedom for the LMS method')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    generate_csv_files()
    estimate_model_all_biomarkers(args)


def generate_csv_files():
    data_file = os.path.join(adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv')
    measurements = vgam.get_measurements_as_collection(data_file)
    for biomarker in adni.biomarker_names:
        print 'Generating output CSV for', biomarker

        csv_file = os.path.join(adni.project_folder, 'data', biomarker.replace(' ', '_') + '.csv')
        writer = csv.writer(open(csv_file, 'wb'), delimiter=',')
        writer.writerow(['rid', 'progress', 'value'])

        for rid, rid_data in measurements.items():
            for _, scan_data in rid_data.items():
                progress = adni.safe_cast(scan_data['progress'], int)
                value = adni.safe_cast(scan_data[biomarker], int)
                if progress != None and value != None:
                    writer.writerow([rid, progress, value])


def estimate_model_all_biomarkers(args):
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(estimate_model)(args, biomarker) for biomarker in adni.biomarker_names)


def estimate_model(args, biomarker):
    print 'Fitting curve to', biomarker
    r_file = os.path.join(adni.project_folder, 'src/fitting/vgam_estimate_curves.R')
    csv_file = os.path.join(adni.project_folder, 'data', biomarker.replace(' ', '_') + '.csv')
    output_file = csv_file.replace('.csv', '_curves.csv')
    densities_file = csv_file.replace('.csv', '_densities.csv')
    image_file = csv_file.replace('.csv', '_curves.pdf')
    stdout_file = csv_file.replace('.csv', '_stdout.Rout')

    command = "R CMD BATCH \"--args input_file='%s' output_file='%s' degrees_of_freedom=%i save_plot=1 plot_file='%s' densities_file='%s'\" %s %s" \
              % (csv_file, output_file, int(args.degrees_of_freedom), image_file, densities_file, r_file, stdout_file)

    subprocess.call(command, shell=True)

if __name__ == '__main__':
    main()
