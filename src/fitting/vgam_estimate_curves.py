#! /usr/bin/env python2.7
import os.path
import argparse
import csv
import subprocess
import joblib as jl
from common import log as log
from common import adni_tools as adni
from vgam.datahandler import DataHandler
from vgam.progressionmodel import MultiBiomarkerProgressionModel
from vgam.modelfitter import ModelFitter


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('-n', '--nr_threads', type=int, default=1, help='number of threads')
    parser.add_argument('-d', '--degrees_of_freedom', type=int, default=2, help='degrees of freedom for the LMS method')
    parser.add_argument('--no_regression', action='store_true', default=False, help='do not perform age regression of biomarker values')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    args = parser.parse_args()

    # Get the data files and biomarkers
    data_handler = DataHandler.get_data_handler(args)

    # Estimate curves
    generate_csv_files(args, data_handler)
    estimate_model_all_biomarkers(args, data_handler)


def generate_csv_files(args, data_handler):
    """
    Generate the CSV file used to call the R script.

    :param Namespace args:
    :param DataHandler data_handler:
    """
    assert isinstance(args, argparse.Namespace)
    assert isinstance(data_handler, DataHandler)

    biomarkers = data_handler.get_biomarker_set()
    measurements = data_handler.get_measurements_as_dict(select_training_set=True)

    if args.iteration > 0:
        scaling_biomarkers = ['CDRSB']
        measurements = data_handler.update_measurements_with_biomarker_values(measurements, biomarkers=scaling_biomarkers)
        model = MultiBiomarkerProgressionModel()
        for biomarker in scaling_biomarkers:
            model_file = data_handler.get_model_file(biomarker, iteration=args.iteration - 1)
            model.add_model(biomarker, model_file)
        measurements = ModelFitter(model).get_scaled_measurements(measurements)

    for biomarker in biomarkers:
        print log.INFO, 'Generating output CSV for {0}...'.format(biomarker)
        samples_file = data_handler.get_samples_file(biomarker)
        writer = csv.writer(open(samples_file, 'wb'), delimiter=',')
        writer.writerow(['rid', 'progress', 'value'])

        subjects = set()
        num_samples = 0
        for rid, visits in measurements.items():
            for _, visit_data in visits.items():
                try:
                    progression = adni.safe_cast(visit_data['progress'], int)
                    biomarker_name = biomarker
                    value = adni.safe_cast(visit_data[biomarker_name], float)
                    if progression is not None and value is not None:
                        writer.writerow([rid, progression, value])
                        subjects.add(rid)
                        num_samples += 1
                except KeyError:
                    pass

        print log.RESULT, 'Collected {0} samples from {1} subjects.'.format(num_samples, len(subjects))


def estimate_model_all_biomarkers(args, data_handler):
    """
    Estimate all models using VGAM.

    :param Namespace args: the command line arguments
    :param DataHandler data_handler: the data handler
    """
    assert isinstance(args, argparse.Namespace)
    assert isinstance(data_handler, DataHandler)

    biomarkers = data_handler.get_biomarker_set()
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(estimate_model)(args, data_handler, biomarker) for biomarker in biomarkers)


def estimate_model(args, data_handler, biomarker):
    """
    Estimate a model using VGAM.

    :param Namespace args: the command line arguments
    :param DataHandler data_handler: the data handler
    :param list biomarker: a list of biomarkers
    """
    assert isinstance(args, argparse.Namespace)
    assert isinstance(data_handler, DataHandler)

    r_file = os.path.join(adni.project_folder, 'src/fitting/vgam_estimate_curves.R')
    samples_file = data_handler.get_samples_file(biomarker)
    model_file = data_handler.get_model_file(biomarker)

    if os.path.isfile(model_file) and not args.recompute_models:
        print log.SKIP, 'Model for {0} already exists.'.format(biomarker)
    else:
        print log.INFO, 'Fitting curve to {0}...'.format(biomarker)
        r_file = os.path.join(adni.project_folder, 'src/fitting/vgam_estimate_curves.R')

        samples_file = data_handler.get_samples_file(biomarker)
        stdout_file = samples_file.replace('_samples.csv', '_stdout.Rout')

        command = "R CMD BATCH \"--args input_file='%s' output_file='%s' degrees_of_freedom=%i \" %s %s" \
                  % (samples_file, model_file, int(args.degrees_of_freedom), r_file, stdout_file)

        subprocess.call(command, shell=True)


if __name__ == '__main__':
    main()
