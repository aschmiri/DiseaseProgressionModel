#! /usr/bin/env python2.7
import os.path
from subprocess import call
import pickle
import numpy as np
from common import log as log
from common.datahandler import SynthDataHandler
from common.progressionmodel import ProgressionModel
from common.modelfitter import ModelFitter
from common.synthmodel import SynthModel


# ------------------------------------------------------------------------------
#  Synthetic data
# ------------------------------------------------------------------------------
def generate_synth_model(biomarker, recompute_models=False, num_samples=1000,
                         sampling='longitudinal', rate_sigma=0.0, conversion_sigma=0.0, run=0):
    # Get model and samples file
    model_file_experiment = SynthDataHandler.get_model_file(biomarker, num_samples=num_samples,
                                                            sampling=sampling, rate_sigma=rate_sigma,
                                                            conversion_sigma=conversion_sigma, run=run)
    samples_file_experiment = SynthDataHandler.get_samples_file(biomarker, num_samples=num_samples,
                                                                sampling=sampling, rate_sigma=rate_sigma,
                                                                conversion_sigma=conversion_sigma, run=run)

    # Read model or recompute
    if os.path.isfile(model_file_experiment) and os.path.isfile(samples_file_experiment) and not recompute_models:
        print log.SKIP, 'Skipping model generation for {0} samples {1}, run {2}'.format(num_samples, sampling, run)
    else:
        src_folder = os.path.join(SynthDataHandler.get_project_folder(), 'src')

        if os.path.isfile(model_file_experiment):
            call(['rm', model_file_experiment])

        while not os.path.isfile(model_file_experiment):
            biomarker_str = ' -b {0}'.format(biomarker) if biomarker is not None else ''
            num_samples_str = ' -n {0}'.format(num_samples) if num_samples is not None else ''
            sampling_str = ' --sampling {0}'.format(sampling)
            rate_sigma_str = ' --rate_sigma {0}'.format(
                rate_sigma) if rate_sigma is not None and rate_sigma > 0.0 else ''
            conv_sigma_str = ' --conversion_sigma {0}'.format(
                conversion_sigma) if conversion_sigma is not None and conversion_sigma > 0.0 else ''

            call('{0}/synth/generate_data.py {1}{2}{3}{4}{5}'.format(
                src_folder, biomarker_str, num_samples_str, sampling_str,
                rate_sigma_str, conv_sigma_str), shell=True)
            call('{0}/training/train_models.py -m synth {1}'.format(
                src_folder, biomarker_str), shell=True)

            model_file = SynthDataHandler.get_model_file(biomarker)
            samples_file = SynthDataHandler.get_samples_file(biomarker)
            if os.path.isfile(model_file) and os.path.isfile(samples_file):
                call(['mv', model_file, model_file_experiment])
                call(['mv', samples_file, samples_file_experiment])
            else:
                print log.WARNING, 'Failed to generate model, retrying...'

    # Return model file
    return model_file_experiment


def evaluate_synth_model(model_file, biomarker, progress_linspace, number_of_value_steps, metric='area'):
    # Define progress steps
    pm = ProgressionModel(biomarker, model_file)
    progresses = np.linspace(progress_linspace[0],
                             progress_linspace[1],
                             progress_linspace[2])

    # Define value steps
    min_val = float('inf')
    max_val = float('-inf')
    for quantile in [0.01, 0.99]:
        curve = pm.get_quantile_curve(progresses, quantile)
        min_val = min(min_val, np.min(curve))
        max_val = max(max_val, np.max(curve))
    values = np.linspace(min_val, max_val, number_of_value_steps)

    # Get mean error
    error = 0
    if metric == 'area':
        for progr in progresses:
            probs_model = [SynthModel.get_probability(biomarker, progr, v) for v in values]
            probs_fit = pm.get_density_distribution(values, progr)
            error += np.sum(np.abs(np.array(probs_fit) - np.array(probs_model)))
        error *= (values[1] - values[0]) / len(progresses)
    elif metric == 'peakdist':
        for progr in progresses:
            probs_model = [SynthModel.get_probability(biomarker, progr, v) for v in values]
            probs_fit = pm.get_density_distribution(values, progr)
            peak_model = values[np.argsort(probs_model)[-1]]
            peak_fit = values[np.argsort(probs_fit)[-1]]
            error += np.abs(peak_fit - peak_model)
        error /= len(progresses)
    elif metric == 'maxdist':
        for value in values:
            probs_model = [SynthModel.get_probability(biomarker, p, value) for p in progresses]
            probs_fit = [pm.get_density_distribution([value], p) for p in progresses]
            max_model = progresses[np.argsort(probs_model)[-1]]
            max_fit = progresses[np.argsort(probs_fit)[-1]]
            error += np.abs(max_fit - max_model)
        error /= len(values)
    else:
        print log.ERROR, 'Metric unknown: {0}'.format(metric)

    return error


def generate_synth_test_data(biomarkers, num_test_samples, number_of_visits, run, recompute_test_data=False):
    biomarkers_str = '_'.join(biomarkers)
    test_data_filename = 'test_data_{0}_{1}_{2}_{3}.p'.format(biomarkers_str, num_test_samples, number_of_visits, run)
    test_data_folder = SynthDataHandler.make_dir(SynthDataHandler.get_eval_folder(), biomarkers_str)
    test_data_file = os.path.join(test_data_folder, test_data_filename)

    if os.path.isfile(test_data_file) and not recompute_test_data:
        print log.SKIP, 'Skipping test data generation with {0} samples...'.format(num_test_samples)
        test_data = pickle.load(open(test_data_file, 'rb'))
    else:
        print log.INFO, 'Generating test data with {0} samples...'.format(num_test_samples)

        test_data = {}
        for rid in xrange(num_test_samples):
            test_data.update({rid: {}})

            base_progress = SynthModel.get_random_progress(sampling='uniform')
            for viscode in xrange(number_of_visits):
                test_data[rid].update({viscode: {}})

                scantime = viscode * 180
                test_data[rid][viscode].update({'scantime': scantime})

                visit_progress = base_progress + scantime
                test_data[rid][viscode].update({'progress': visit_progress})

                for biomarker in biomarkers:
                    value = SynthModel.get_distributed_value(biomarker, visit_progress)
                    test_data[rid][viscode].update({biomarker: value})
        pickle.dump(test_data, open(test_data_file, 'wb'))

    return test_data


def evaluate_synth_fitting(fitter, test_data, biomarkers, viscodes):
    assert isinstance(fitter, ModelFitter)

    print log.INFO, 'Testing biomarkers {0} with viscodes {1}...'.format(biomarkers, viscodes)
    dpis = []
    progresses = []

    for rid in test_data:
        # Get progress
        progress = test_data[rid][0]['progress']

        # Estimate DPI
        samples = {}
        for viscode in viscodes:
            samples.update({viscode: {}})
            samples[viscode].update({'scantime': test_data[rid][viscode]['scantime']})
            for biomarker in biomarkers:
                samples[viscode].update({biomarker: test_data[rid][viscode][biomarker]})
        dpi = fitter.get_dpi_for_samples(samples)

        if dpi is not None:
            dpis.append(dpi)
            progresses.append(progress)

    errors = np.abs(np.array(progresses) - np.array(dpis))
    rms_error = np.sqrt(np.sum(np.square(errors)) / len(dpis))
    mean_error = np.sum(errors) / len(dpis)
    print log.RESULT, 'Mean error: {0}, RMS error: {1}'.format(mean_error, rms_error)

    return errors
