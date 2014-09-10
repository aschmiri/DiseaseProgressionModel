#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
from common import log as log
from common import adni_tools as adni
from vgam.datahandler import SynthDataHandler
from vgam.synthmodel import SynthModel
from vgam.progressionmodel import MultiBiomarkerProgressionModel
from vgam.modelfitter import ModelFitter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--recompute_errors', action='store_true', help='recompute the errors of the models')
    parser.add_argument('-n', '--number_of_samples', type=int, default=100, help='the number of test samples')
    args = parser.parse_args()

    # Initialise model fitter
    data_handler = SynthDataHandler()
    biomarkers = ['synth_e1', 'synth_e2']  # SynthModel.get_biomarker_names()

    model = MultiBiomarkerProgressionModel()
    for biomarker in biomarkers:
        model_file = data_handler.get_model_file(biomarker)
        model.add_model(biomarker, model_file)
    fitter = ModelFitter(model)

    if args.recompute_errors:
        # Generate test samples
        test_data = generate_test_data(biomarkers, args.number_of_samples, 5)

        # Test different biomarker constellations
        errors = []
#         errors.append(test_fitting(test_data, fitter, [0], ['synth0']))
#         errors.append(test_fitting(test_data, fitter, [0], ['synth1']))
#         errors.append(test_fitting(test_data, fitter, [0], ['synth2']))
#         errors.append(test_fitting(test_data, fitter, [0], ['synth0', 'synth1', 'synth2']))
#
#         errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth0']))
#         errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth1']))
#         errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth2']))
#         errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth0', 'synth1', 'synth2']))
#
#         errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth0']))
#         errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth1']))
#         errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth2']))
#         errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth0', 'synth1', 'synth2']))

        errors.append(test_fitting(test_data, fitter, [0], ['synth_e1']))
        errors.append(test_fitting(test_data, fitter, [0], ['synth_e2']))
        errors.append(test_fitting(test_data, fitter, [0], ['synth_e1', 'synth_e2']))

        errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth_e1']))
        errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth_e2']))
        errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth_e1', 'synth_e2']))

        errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth_e1']))
        errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth_e2']))
        errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth_e1', 'synth_e2']))

        save_errors(errors)
    else:
        errors = read_errors()

    plot_errors(errors)


def generate_test_data(biomarkers, number_of_samples, number_of_visits):
    print log.INFO, 'Generating test set with {0} samples...'.format(number_of_samples)

    test_data = {}
    for rid in xrange(number_of_samples):
        test_data.update({rid: {}})

        base_progress = SynthModel.get_random_progress(uniform_progression=True)
        for viscode in xrange(number_of_visits):
            test_data[rid].update({viscode: {}})

            scantime = viscode * 180
            test_data[rid][viscode].update({'scantime': scantime})

            visit_progress = base_progress + scantime
            test_data[rid][viscode].update({'progress': visit_progress})

            for biomarker in biomarkers:
                value = SynthModel.get_distributed_value(biomarker, visit_progress)
                test_data[rid][viscode].update({biomarker: value})

    return test_data


def test_fitting(test_data, fitter, test_viscodes, test_biomarkers):
    print log.INFO, 'Testing biomarkers {0} with viscodes {1}...'.format(test_biomarkers, test_viscodes)
    dpis = []
    progresses = []

    for rid in test_data:
        # Get progress
        progress = test_data[rid][0]['progress']

        # Estimate DPI
        samples = {}
        for viscode in test_viscodes:
            samples.update({viscode: {}})
            samples[viscode].update({'scantime': test_data[rid][viscode]['scantime']})
            for biomarker in test_biomarkers:
                samples[viscode].update({biomarker: test_data[rid][viscode][biomarker]})
        dpi = fitter.get_dpi_for_samples(samples)

        # print log.RESULT, 'Estimated DPI: {0}, Progress: {1}'.format(dpi, progress)
        if dpi is not None:
            dpis.append(dpi)
            progresses.append(progress)

    errors = np.abs(np.array(progresses) - np.array(dpis))
    rms_error = np.sqrt(np.sum(np.square(errors)) / len(dpis))
    mean_error = np.sum(errors) / len(dpis)
    print log.RESULT, 'Mean error: {0}, RMS error: {1}'.format(mean_error, rms_error)

    return errors


def read_errors():
    evaluation_file = os.path.join(adni.eval_folder, 'eval_synth_fitting_e.p')
    return pickle.load(open(evaluation_file, 'rb'))


def save_errors(errors):
    evaluation_file = os.path.join(adni.eval_folder, 'eval_synth_fitting_e.p')
    pickle.dump(errors, open(evaluation_file, 'wb'))


def plot_errors(errors):
    plt.boxplot(errors)
    plt.show()

if __name__ == '__main__':
    main()
