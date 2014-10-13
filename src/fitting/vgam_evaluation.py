#! /usr/bin/env python2.7
import os.path
from subprocess import call
import pickle
import numpy as np
import matplotlib as mpl
from common import log as log
from common import adni_tools as adni
from vgam.synthmodel import SynthModel
from vgam.datahandler import DataHandler
from vgam.progressionmodel import ProgressionModel
from vgam.progressionmodel import MultiBiomarkerProgressionModel
from vgam.modelfitter import ModelFitter


def get_estimates_file(method, visits, biomarkers, estimate_dprs):
    estimates_file_trunk = 'estimate_dpi_dpr_with_{0}_{1}.p' if estimate_dprs else 'estimate_dpi_with_{0}_{1}.p'
    if biomarkers is None:
        estimates_file_basename = estimates_file_trunk.format(method, '_'.join(visits))
    else:
        biomarkers_string = '_'.join(biomarkers).replace(' ', '_')
        estimates_file_basename = estimates_file_trunk.format(biomarkers_string, '_'.join(visits))
    estimates_file = os.path.join(adni.eval_folder, estimates_file_basename)
    return estimates_file


def get_prediction_file(method, visits, biomarkers, estimate_dprs, use_last_visit, predict_biomarker):
    predict_biomarker_str = predict_biomarker.replace(' ', '_')
    predict_file_trunk = 'predict_{0}_with_dpr_{1}_{2}{3}.p' if estimate_dprs else 'predict_{0}_with_{1}_{2}{3}.p'
    if biomarkers is None:
        predict_file_basename = predict_file_trunk.format(predict_biomarker_str,
                                                          method, '_'.join(visits),
                                                          '_last' if use_last_visit else '')
    else:
        estimate_biomarkers_string = '_'.join(biomarkers).replace(' ', '_')
        predict_file_basename = predict_file_trunk.format(predict_biomarker_str,
                                                          estimate_biomarkers_string,
                                                          '_'.join(visits),
                                                          '_last' if use_last_visit else '')
    prediction_file = os.path.join(adni.eval_folder, predict_file_basename)

    return prediction_file


def get_progress_estimates(args, estimates_file=None):
    # Get filename
    if estimates_file is None:
        estimates_file = get_estimates_file(args.method, args.visits, args.biomarkers, args.estimate_dprs)

    # Read if estimates exist, else recompute
    if os.path.isfile(estimates_file) and not args.recompute_estimates:
        # Read test results from file
        print log.INFO, 'Reading DPI{0} estimations from {1}...'.format('\DPR' if args.estimate_dprs else '', estimates_file)
        (rids, diagnoses, dpis, dprs, mean_min, mean_max) = pickle.load(open(estimates_file, 'rb'))
    else:
        # Collect data for test
        data_handler = DataHandler.get_data_handler(args)
        biomarkers = data_handler.get_biomarker_set()
        measurements = data_handler.get_measurements_as_dict(visits=args.visits,
                                                             biomarkers=biomarkers,
                                                             select_test_set=True,
                                                             select_complete=True)

        # Setup model
        model = MultiBiomarkerProgressionModel()
        for biomarker in biomarkers:
            model_file = data_handler.get_model_file(biomarker)
            model.add_model(biomarker, model_file)
        fitter = ModelFitter(model)

        # Calculate mean and max progress
        mean_min = model.get_mean_min_progress()
        mean_max = model.get_mean_max_progress()

        # Estimate dpis (and dprs) and save data
        if not args.estimate_dprs or len(args.visits) == 1:
            if args.estimate_dprs and len(args.visits) == 1:
                print log.WARNING, 'Only one visit, cannot estimate DPR (setting to one)'
            rids, diagnoses, dpis = estimate_dpis(measurements, args.visits, fitter)
            dprs = np.ones(len(dpis)).tolist()
        else:
            rids, diagnoses, dpis, dprs = estimate_dpis_dprs(measurements, args.visits, fitter)

        print log.INFO, 'Saving DPI{0} estimations to {1}...'.format('\DPR' if args.estimate_dprs else '', estimates_file)
        pickle.dump((rids, diagnoses, dpis, dprs, mean_min, mean_max), open(estimates_file, 'wb'))

    # Reduce to consistent data sets with bl, m12 and m24 samples
    if args.consistent_data:
        data_handler = DataHandler.get_data_handler(args)
        consistent_measures = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24'],
                                                                    select_test_set=True,
                                                                    select_complete=True)
        consistent_rids = []
        consistent_diagnoses = []
        consistent_dpis = []
        consistent_dprs = []
        for i, rid in enumerate(rids):
            if rid in consistent_measures:
                consistent_rids.append(rid)
                consistent_diagnoses.append(diagnoses[i])
                consistent_dpis.append(dpis[i])
                consistent_dprs.append(dprs[i])
        rids = consistent_rids
        diagnoses = consistent_diagnoses
        dpis = consistent_dpis
        dprs = consistent_dprs

        print log.INFO, 'Selected {0} consistent subjects.'.format(len(dpis))

    # Return results
    return rids, diagnoses, dpis, dprs, mean_min, mean_max


def estimate_dpis(measurements, viscodes, fitter):
    rids = []
    diagnoses = []
    dpis = []
    for rid in measurements:
        print log.INFO, 'Estimating DPI for subject {0}...'.format(rid)
        try:
            diagnosis = measurements[rid]['bl']['DX.scan']
        except KeyError:
            continue

        if not set(viscodes).issubset(set(measurements[rid].keys())):
            print log.WARNING, 'Not all viscodes {0} available for subject {1}!'.format(viscodes, rid)
            continue

        samples = {}
        for viscode in viscodes:
            samples.update({viscode: measurements[rid][viscode]})
        dpi = fitter.get_dpi_for_samples(samples)

        print log.RESULT, 'Subject {0}: Estimated DPI: {1}, Diagnosis: {2}'.format(rid, dpi, diagnosis)
        if dpi is not None:
            rids.append(rid)
            diagnoses.append(diagnosis)
            dpis.append(dpi)

    return rids, diagnoses, dpis


def estimate_dpis_dprs(measurements, viscodes, fitter):
    # Test all available subjects
    rids = []
    diagnoses = []
    dpis = []
    dprs = []
    for rid in measurements:
        print log.INFO, 'Estimating DPI and DPR for subject {0}...'.format(rid)
        try:
            diagnosis = measurements[rid]['bl']['DX.scan']
        except KeyError:
            continue

        if not set(viscodes).issubset(set(measurements[rid].keys())):
            print log.WARNING, 'Not all viscodes {0} available for subject {1}!'.format(viscodes, rid)
            continue

        samples = {}
        for viscode in viscodes:
            samples.update({viscode: measurements[rid][viscode]})
        dpi, dpr = fitter.get_dpi_dpr_for_samples(samples)

        print log.RESULT, 'Subject {0}: Estimated DPI: {1}, DPR: {2}, Diagnosis: {3}'.format(rid, dpi, dpr, diagnosis)
        if dpi is not None and dpr is not None:
            rids.append(rid)
            diagnoses.append(diagnosis)
            dpis.append(dpi)
            dprs.append(dpr)

    # Plot the results
    return rids, diagnoses, dpis, dprs


def generate_synth_model(args, data_handler, biomarker, num_samples=1000, sampling='longitudinal', rate_sigma=0.0, run=0):
    # Get model and samples file
    model_file_experiment = data_handler.get_model_file(biomarker, num_samples=num_samples, sampling=sampling,
                                                        rate_sigma=rate_sigma, run=run)
    samples_file_experiment = data_handler.get_samples_file(biomarker, num_samples=num_samples, sampling=sampling,
                                                            rate_sigma=rate_sigma, run=run)

    # Read model or recompute
    if os.path.isfile(model_file_experiment) and os.path.isfile(samples_file_experiment) and not args.recompute_models:
        print log.SKIP, 'Skipping model generation for {0} samples {1}, run {2}'.format(num_samples, sampling, run)
    else:
        exec_folder = os.path.join(adni.project_folder, 'src', 'fitting')

        if os.path.isfile(model_file_experiment):
            call(['rm', model_file_experiment])

        while not os.path.isfile(model_file_experiment):
            biomarker_str = ' -b {0}'.format(biomarker) if biomarker is not None else ''
            num_samples_str = ' -n {0}'.format(num_samples) if num_samples is not None else ''
            sampling_str = ' --sampling {0}'.format(sampling)
            rate_sigma_str = ' --rate_sigma {0}'.format(
                rate_sigma) if rate_sigma is not None and rate_sigma > 0.0 else ''

            call('{0}/vgam_generate_synth_data.py {1}{2}{3}{4}'.format(exec_folder, biomarker_str, num_samples_str,
                                                                       sampling_str, rate_sigma_str), shell=True)
            call('{0}/vgam_estimate_curves.py synth {1}'.format(exec_folder, biomarker_str), shell=True)

            model_file = data_handler.get_model_file(biomarker)
            samples_file = data_handler.get_samples_file(biomarker)
            if os.path.isfile(model_file) and os.path.isfile(samples_file):
                call(['mv', model_file, model_file_experiment])
                call(['mv', samples_file, samples_file_experiment])
            else:
                print log.WARNING, 'Failed to generate model, retrying...'

    # Return model file
    return model_file_experiment


def evaluate_synth_model(args, model_file, biomarker, metric='area'):
    # Define progress steps
    pm = ProgressionModel(biomarker, model_file)
    progresses = np.linspace(args.progress_linspace[0],
                             args.progress_linspace[1],
                             args.progress_linspace[2])

    # Define value steps
    min_val = float('inf')
    max_val = float('-inf')
    for quantile in [0.01, 0.99]:
        curve = pm.get_quantile_curve(progresses, quantile)
        min_val = min(min_val, np.min(curve))
        max_val = max(max_val, np.max(curve))
    values = np.linspace(min_val, max_val, args.number_of_value_steps)

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


def generate_synth_test_data(args, biomarkers, num_test_samples, number_of_visits, run=None):
    print log.INFO, 'Generating test set with {0} samples...'.format(num_test_samples)
    biomarkers_str = '_'.join(biomarkers)
    if run is None:
        test_data_filename = 'test_data_{0}_{1}_{2}.p'.format(biomarkers_str, num_test_samples, number_of_visits)
    else:
        test_data_filename = 'test_data_{0}_{1}_{2}_{3}.p'.format(biomarkers_str, num_test_samples, number_of_visits,
                                                                  run)
    test_data_file = os.path.join(adni.eval_folder, test_data_filename)

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


def evaluate_synth_fitting(fitter, test_data, biomarkers, viscodes=None):
    assert isinstance(fitter, ModelFitter)
    viscodes = [0] if viscodes is None else viscodes

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
    assert isinstance(ax, mpl.axes.Axes)

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.xaxis.grid(True, linestyle=':', which='major', color='lightgrey', alpha=0.5)
    ax.yaxis.grid(True, linestyle=':', which='major', color='lightgrey', alpha=0.5)
    ax.tick_params(axis='both', which='both', bottom='on', top='off', left='off', right='off')

    plt.setp(ax.get_yticklabels(), rotation=45, fontsize=10)


def get_metric_unit(biomarker):
    if biomarker in adni.structure_names + ['synth_hipp', 'synth_brain']:
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
