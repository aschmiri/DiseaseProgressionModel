#! /usr/bin/env python2.7
import os.path
import argparse
import csv
import joblib as jl
from subprocess import call
from common import log as log
from common.datahandler import DataHandler


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-n', '--nr_threads', type=int, default=1, help='number of threads')
    parser.add_argument('-d', '--degrees_of_freedom', type=int, default=2, help='degrees of freedom for the LMS method')
    parser.add_argument('--min_visits', type=int, default=0, help='the minimal number of visits')
    parser.add_argument('--no_regression', action='store_true', default=False, help='do not perform age regression of biomarker values')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    args = parser.parse_args()

    # Get the data files and biomarkers
    data_handler = DataHandler.get_data_handler(method=args.method, biomarkers=args.biomarkers)

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

    biomarkers = data_handler.get_biomarker_names()
    measurements = data_handler.get_measurements_as_dict(min_visits=args.min_visits, select_training_set=True)

    for biomarker in biomarkers:
        print log.INFO, 'Generating output CSV for {0}...'.format(biomarker)
        samples_file = data_handler.get_samples_file(biomarker)
        writer = csv.writer(open(samples_file, 'wb'), delimiter=',')
        writer.writerow(['rid', 'progress', 'value', 'diagnosis'])

        subjects = set()
        num_samples = 0
        for rid, visits in measurements.items():
            for _, visit_data in visits.items():
                try:
                    progress = DataHandler.safe_cast(visit_data['progress'], int)
                    biomarker_name = biomarker
                    value = DataHandler.safe_cast(visit_data[biomarker_name], float)
                    diagnosis = DataHandler.safe_cast(visit_data['DX.scan'], float)
                    if progress is not None and value is not None:
                        writer.writerow([rid, progress, value, diagnosis])
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

    biomarkers = data_handler.get_biomarker_names()
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(estimate_model)(args, biomarker) for biomarker in biomarkers)


def estimate_model(args, biomarker):
    """
    Estimate a model using VGAM.

    :param Namespace args: the command line arguments
    :param list biomarker: a list of biomarkers
    """
    assert isinstance(args, argparse.Namespace)

    model_file = DataHandler.get_model_file(biomarker)
    if os.path.isfile(model_file) and not args.recompute_models:
        print log.SKIP, 'Model for {0} already exists.'.format(biomarker)
    else:
        # Clean old model
        if os.path.isfile(model_file):
            call(['rm', model_file])

        # Estiamte model
        print log.INFO, 'Fitting curve to {0}...'.format(biomarker)
        r_file = os.path.join(DataHandler.get_project_folder(), 'src/training/train_models.R')
        samples_file = DataHandler.get_samples_file(biomarker)
        stdout_file = samples_file.replace('_samples.csv', '_stdout.Rout')

        command = "R CMD BATCH \"--args input_file='%s' output_file='%s' degrees_of_freedom=%i \" %s %s" \
                  % (samples_file, model_file, int(args.degrees_of_freedom), r_file, stdout_file)

        call(command, shell=True)

        # Check if model was generated
        if not os.path.isfile(model_file):
            print log.ERROR, 'Failed to generate model for {0}!'.format(biomarker)


if __name__ == '__main__':
    main()