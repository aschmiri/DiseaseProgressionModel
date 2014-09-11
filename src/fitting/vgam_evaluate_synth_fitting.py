#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
from subprocess import call
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
    parser.add_argument('-n', '--number_of_test_samples', type=int, default=100, help='the number of test samples')
    parser.add_argument('--recompute_errors', action='store_true', help='recompute the errors of the models')
    parser.add_argument('--recompute_models', action='store_true', help='recompute the models with new samples')
    parser.add_argument('--number_of_training_samples', type=int, default=1000, help='the number of training samples')
    parser.add_argument('--output_file', type=str, default=None, help='filename of the output image with the plot')
    args = parser.parse_args()

    evaluation_file = os.path.join(adni.eval_folder, 'eval_synth_fitting_e.p')
    if os.path.isfile(evaluation_file) and not args.recompute_errors:
        # Read test results from file
        print log.INFO, 'Reading test results from {0}...'.format(evaluation_file)
        errors = pickle.load(open(evaluation_file, 'rb'))
    else:
        # Generate test results
        print log.INFO, 'Starting tests...'

        # Initialise data handler
        data_handler = SynthDataHandler()
        biomarkers = SynthModel.get_biomarker_names()

        # Initialise model fitter
        model = MultiBiomarkerProgressionModel()
        for biomarker in biomarkers:
            model_file = data_handler.get_model_file(biomarker, num_samples=args.number_of_training_samples)
            if os.path.isfile(model_file) and not args.recompute_models:
                print log.SKIP, 'Skipping model generation for {0} with {1} samples'.format(biomarker, args.number_of_training_samples)
            else:
                generate_model(args, data_handler, biomarker, args.number_of_training_samples)
            model.add_model(biomarker, model_file)
        fitter = ModelFitter(model)

        # Generate test data and test different biomarker constellations
        test_data = generate_test_data(biomarkers, args.number_of_test_samples, 5)
        errors = run_test_1(test_data, fitter)

        # Save results
        print log.INFO, 'Saving test results to {0}...'.format(evaluation_file)
        pickle.dump(errors, open(evaluation_file, 'wb'))

    # Analyse and plot the test results
    plot_errors(args, errors)


def generate_model(args, data_handler, biomarker, num_training_samples):
    model_file = data_handler.get_model_file(biomarker)
    model_file_experiment = data_handler.get_model_file(biomarker, num_samples=num_training_samples)
    samples_file = data_handler.get_samples_file(biomarker)
    samples_file_experiment = data_handler.get_samples_file(biomarker, num_samples=num_training_samples)

    exec_folder = os.path.join(adni.project_folder, 'src', 'fitting')
    call('{0}/vgam_generate_synth_data.py -b {1} -n {2}'.format(exec_folder, biomarker, num_training_samples), shell=True)
    call('{0}/vgam_estimate_curves.py synth -b {1}'.format(exec_folder, biomarker), shell=True)
    call(['mv', model_file, model_file_experiment])
    call(['mv', samples_file, samples_file_experiment])


def generate_test_data(biomarkers, num_test_samples, number_of_visits):
    print log.INFO, 'Generating test set with {0} samples...'.format(num_test_samples)

    test_data = {}
    for rid in xrange(num_test_samples):
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


def run_test_1(test_data, fitter):
    errors = []

    errors.append(test_fitting(test_data, fitter, [0], ['synth0']))
    errors.append(test_fitting(test_data, fitter, [0], ['synth1']))
    errors.append(test_fitting(test_data, fitter, [0], ['synth2']))
    errors.append(test_fitting(test_data, fitter, [0], ['synth0', 'synth1', 'synth2']))

    errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth0']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth1']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth2']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth0', 'synth1', 'synth2']))

    errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth0']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth1']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth2']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth0', 'synth1', 'synth2']))

    return errors


def run_test_2(test_data, fitter):
    errors = []

    errors.append(test_fitting(test_data, fitter, [0], ['synth_e1']))
    errors.append(test_fitting(test_data, fitter, [0], ['synth_e2']))
    errors.append(test_fitting(test_data, fitter, [0], ['synth_e1', 'synth_e2']))

    errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth_e1']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth_e2']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2], ['synth_e1', 'synth_e2']))

    errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth_e1']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth_e2']))
    errors.append(test_fitting(test_data, fitter, [0, 1, 2, 3, 4], ['synth_e1', 'synth_e2']))

    return errors


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


def plot_errors(args, errors):
    print log.INFO, 'Plotting results...'

    # Define figure and grid
    fig, ax1 = plt.subplots(figsize=(10, 6))
    top = 1800
    bottom = -50
    ax1.set_ylim(bottom, top)
    ax1.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.7)
    ax1.axvline(3.5, color='lightgrey', alpha=0.7)
    ax1.axvline(6.5, color='lightgrey', alpha=0.7)

    # Set labels
    ax1.set_axisbelow(True)
    ax1.set_title('DPI estimation using different biomarker settings')
    ax1.set_ylabel('Fitting error')
    ax1.set_xticklabels(['S1', 'S2', 'S1 + S2', 'S1', 'S2', 'S1 + S2', 'S1', 'S2', 'S1 + S2'])
    ax1.text(2, -220, '1 timepoint', horizontalalignment='center')
    ax1.text(5, -220, '3 timepoints', horizontalalignment='center')
    ax1.text(8, -220, '5 timepoints', horizontalalignment='center')

    # Write mean errors as text
    pos = np.arange(len(errors)) + 1
    medians = [np.median(e) for e in errors]
    upperLabels = [str(np.round(s, 2)) for s in medians]
    weights = ['semibold', 'semibold', 'bold'] * 3
    for tick in range(len(errors)):
        ax1.text(pos[tick], 0.95 * top, upperLabels[tick],
                 horizontalalignment='center', size='x-small', weight=weights[tick])

    # Plot boxplots
    plt.boxplot(errors)

    # Show or save plot
    if args.output_file is not None:
        plt.savefig(args.output_file)
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
