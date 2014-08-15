#! /usr/bin/env python2.7
import os.path
import argparse
import csv
import subprocess
import joblib as jl
import vgam
from bmia.common import log as log
from bmia.common import adni_tools as adni


def main():
    # parse command line options
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser.add_argument('method', choices=['reg', 'long', 'cons', 'graph'])
    parser.add_argument('-n', '--nr_threads', dest='nr_threads', type=int, default=4, help='number of threads')
    parser.add_argument('-d', '--degrees_of_freedom', dest='degrees_of_freedom', type=int, default=2, help='degrees of freedom for the LMS method')
    parser.add_argument('-s', '--scale_measurements', dest='scale_measurements', action='store_true', default=False, help='scale the measurements by fitting to an initial model')
    parser.add_argument('--trans', dest='trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic (regbased only)')
    parser.add_argument('--spacing', type=str, default='5', help='the transformation spacing (regbased only)')
    args = parser.parse_args()

    # setup data file
    data_filename = {'reg': 'volumes_segbased_' + args.trans + '_' + args.spacing + 'mm.csv',
                     'long': 'volumes_segbased_longitudinal.csv',
                     'cons': 'volumes_segbased_consistent.csv',
                     'graph': 'volumes_segbased_graphcut.csv'}
    data_filepath = os.path.join(adni.project_folder, 'lists', data_filename[args.method])

    # Set the data file and the biomarkers to be considered
    biomarkers = adni.volume_names_essential

    # Estimate curves
    generate_csv_files(args, biomarkers, data_filepath)
    estimate_model_all_biomarkers(args, biomarkers)


def generate_csv_files(args, biomarkers, data_file):
    out_folder = os.path.join(adni.project_folder, 'data', args.method)
    measurements = vgam.get_measurements_as_collection(data_file)
    if args.scale_measurements:
        measurements = vgam.get_scaled_measurements(measurements, biomarkers=['CDRSB'])

    for biomarker in biomarkers:
        print log.INFO, 'Generating output CSV for {0}...'.format(biomarker)

        csv_file = os.path.join(out_folder, biomarker.replace(' ', '_') + '.csv')
        writer = csv.writer(open(csv_file, 'wb'), delimiter=',')
        writer.writerow(['rid', 'progress', 'value'])

        for rid, rid_data in measurements.items():
            for _, scan_data in rid_data.items():
                progress = adni.safe_cast(scan_data['progress'], int)
                value = adni.safe_cast(scan_data[biomarker], int)
                if progress is not None and value is not None and value >= 0:
                    writer.writerow([rid, progress, value])


def estimate_model_all_biomarkers(args, biomarkers):
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(estimate_model)(args, biomarker) for biomarker in biomarkers)


def estimate_model(args, biomarker):
    print log.INFO, 'Fitting curve to {0}...'.format(biomarker)
    r_file = os.path.join(adni.project_folder, 'src/bmia/fitting/vgam_estimate_curves.R')
    out_folder = os.path.join(adni.project_folder, 'data', args.method)
    csv_file = os.path.join(out_folder, biomarker.replace(' ', '_') + '.csv')
    output_file = csv_file.replace('.csv', '_curves.csv')
    densities_file = csv_file.replace('.csv', '_densities.csv')
    image_file = csv_file.replace('.csv', '_curves.pdf')
    stdout_file = csv_file.replace('.csv', '_stdout.Rout')

    command = "R CMD BATCH \"--args input_file='%s' output_file='%s' degrees_of_freedom=%i save_plot=1 plot_file='%s' densities_file='%s'\" %s %s" \
              % (csv_file, output_file, int(args.degrees_of_freedom), image_file, densities_file, r_file, stdout_file)

    subprocess.call(command, shell=True)


if __name__ == '__main__':
    main()
