#! /usr/bin/env python2.7
import os.path
from subprocess import call
import pickle
import numpy as np
from common import log as log
from common import adni_tools as adni
from vgam.synthmodel import SynthModel
from vgam.progressionmodel import ProgressionModel


def generate_model(args, data_handler, biomarker, num_samples=1000, sampling='longitudinal', run=0):
    model_file_experiment = data_handler.get_model_file(biomarker, num_samples=num_samples, sampling=sampling, run=run)
    samples_file_experiment = data_handler.get_samples_file(biomarker, num_samples=num_samples, sampling=sampling, run=run)

    if os.path.isfile(model_file_experiment) and os.path.isfile(samples_file_experiment) and not args.recompute_models:
        print log.SKIP, 'Skipping model generation for {0} samples {1}, run {2}'.format(num_samples, sampling, run)
    else:
        exec_folder = os.path.join(adni.project_folder, 'src', 'fitting')

        if os.path.isfile(model_file_experiment):
            call(['rm', model_file_experiment])

        while not os.path.isfile(model_file_experiment):
            call('{0}/vgam_generate_synth_data.py -b {1} -n {2} --sampling {3}'.format(exec_folder, biomarker, num_samples, sampling), shell=True)
            call('{0}/vgam_estimate_curves.py synth -b {1}'.format(exec_folder, biomarker), shell=True)

            model_file = data_handler.get_model_file(biomarker)
            samples_file = data_handler.get_samples_file(biomarker)
            if os.path.isfile(model_file) and os.path.isfile(samples_file):
                call(['mv', model_file, model_file_experiment])
                call(['mv', samples_file, samples_file_experiment])
            else:
                print log.WARNING, 'Failed to generate model, retrying...'
    return model_file_experiment


def evaluate_synth_model(args, model_file, biomarker):

    # Define progression steps
    pm = ProgressionModel(biomarker, model_file)
    progressions = np.linspace(-args.progression_range,
                               args.progression_range,
                               args.number_of_progression_steps)

    # Define value steps
    min_val = float('inf')
    max_val = float('-inf')
    for std in [-2.0, 2.0]:
        curve = pm.get_quantile_curve(progressions, std)
        min_val = min(min_val, np.min(curve))
        max_val = max(max_val, np.max(curve))
    values = np.linspace(min_val, max_val, args.number_of_value_steps)

    # Get mean error
    error = 0
    for progr in progressions:
        probs_fit = pm.get_density_distribution(values, progr)
        probs_model = [SynthModel.get_probability(biomarker, progr, v) for v in values]
        error += np.sum(np.abs(np.array(probs_fit) - np.array(probs_model)))
    error *= (values[1] - values[0]) / len(progressions)

    return error


def generate_test_data(args, biomarkers, num_test_samples, number_of_visits):
    print log.INFO, 'Generating test set with {0} samples...'.format(num_test_samples)
    biomarkers_str = '_'.join(biomarkers)
    test_data_file = os.path.join(adni.eval_folder, 'test_data_{0}_{1}_{2}.p'.format(biomarkers_str, num_test_samples, number_of_visits))

    if os.path.isfile(test_data_file) and not args.recompute_test_data:
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
    else:
        test_data = pickle.load(open(test_data_file, 'rb'))

    return test_data


def evaluate_synth_fitting(fitter, test_data, biomarkers, viscodes=[0]):
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


def setup_axes(plt, ax):
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.xaxis.grid(True, linestyle=':', which='major', color='lightgrey', alpha=0.5)
    ax.yaxis.grid(True, linestyle=':', which='major', color='lightgrey', alpha=0.5)
    ax.tick_params(axis='both', which='both', bottom='on', top='off', left='off', right='off')

    plt.setp(ax.get_yticklabels(), rotation=45, fontsize=10)


def get_metric_unit(biomarker):
    if biomarker in adni.volume_names + ['synth_hipp', 'synth_brain']:
        return 'Volume'
    if biomarker in adni.cog_score_names + ['synth_mmse', 'synth_cdrsb']:
        return 'Score'
    if biomarker in adni.manifold_coordinate_names:
        return 'Value'
    else:
        return None


def set_boxplot_color(boxplot, index, color):
    box = boxplot['boxes'][index]
    box.set_facecolor(color + (0.2,))
    box.set_edgecolor(color)
    median = boxplot['medians'][index]
    median.set_color(color)
    cap = boxplot['caps'][2 * index]
    cap.set_color(color)
    cap = boxplot['caps'][2 * index + 1]
    cap.set_color(color)
    whisker = boxplot['whiskers'][2 * index]
    whisker.set_linestyle('-')
    whisker.set_color(color)
    whisker = boxplot['whiskers'][2 * index + 1]
    whisker.set_linestyle('-')
    whisker.set_color(color)
    flier = boxplot['fliers'][2 * index]
    flier.set_color(color)
    flier = boxplot['fliers'][2 * index + 1]
    flier.set_color(color)
