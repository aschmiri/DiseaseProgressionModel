#! /usr/bin/env python2.7
import os.path
from subprocess import call
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from common import log as log
from common import adni_tools as adni
from common import adni_plot as aplt
from vgam.synthmodel import SynthModel
from vgam.datahandler import DataHandler
from vgam.progressionmodel import ProgressionModel
from vgam.progressionmodel import MultiBiomarkerProgressionModel
from vgam.modelfitter import ModelFitter


# ------------------------------------------------------------------------------
#  Progression estimation
# ------------------------------------------------------------------------------
def get_progress_estimates(args, visits, method=None, recompute_estimates=False,
                           estimate_dprs=False, consistent_data=False):
    # Overwrite method if defined
    method = method if method is not None else args.method

    # Get filename
    estimates_file_trunk = 'estimate_dpi_dpr_with_{0}_{1}.p' if estimate_dprs else 'estimate_dpi_with_{0}_{1}.p'
    if args.biomarkers is None:
        estimates_file_basename = estimates_file_trunk.format(method, '_'.join(visits))
    else:
        biomarkers_string = '_'.join(args.biomarkers).replace(' ', '_')
        estimates_file_basename = estimates_file_trunk.format(biomarkers_string, '_'.join(visits))
    estimates_file = os.path.join(adni.eval_folder, estimates_file_basename)

    # Read if estimates exist, else recompute
    if os.path.isfile(estimates_file) and not recompute_estimates:
        # Read test results from file
        print log.INFO, 'Reading DPI{0} estimations from {1}...'.format('\DPR' if estimate_dprs else '', estimates_file)
        (rids, diagnoses, dpis, dprs, mean_min, mean_max) = pickle.load(open(estimates_file, 'rb'))
    else:
        # Collect data for test
        data_handler = DataHandler.get_data_handler(args)
        biomarkers = data_handler.get_biomarker_set()
        measurements = data_handler.get_measurements_as_dict(visits=visits,
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
        if not estimate_dprs or len(visits) == 1:
            if estimate_dprs and len(visits) == 1:
                print log.WARNING, 'Only one visit, cannot estimate DPR (setting to one)'
            rids, diagnoses, dpis = estimate_dpis(measurements, visits, fitter)
            dprs = np.ones(len(dpis)).tolist()
        else:
            rids, diagnoses, dpis, dprs = estimate_dpis_dprs(measurements, visits, fitter)

        print log.INFO, 'Saving DPI{0} estimations to {1}...'.format('\DPR' if estimate_dprs else '', estimates_file)
        pickle.dump((rids, diagnoses, dpis, dprs, mean_min, mean_max), open(estimates_file, 'wb'))

    # Reduce to consistent data sets with bl, m12 and m24 samples
    if consistent_data:
        data_handler = DataHandler.get_data_handler(args)
        consistent_measures = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24'],
                                                                    biomarkers=adni.biomarker_names,
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


# ------------------------------------------------------------------------------
#  Biomarker predictions
# ------------------------------------------------------------------------------
def get_biomarker_predictions(args, visits, biomarker, method=None,
                              recompute_estimates=False, recompute_predictions=False,
                              estimate_dprs=False, consistent_data=False, exclude_cn=False,
                              use_last_visit=False, naive_use_diagnosis=False,
                              plot_predictions=False):
    method = method if method is not None else args.method

    # Get prediction file
    predict_biomarker_str = biomarker.replace(' ', '_')
    predict_file_trunk = 'predict_{0}_with_dpr_{1}_{2}{3}.p' if estimate_dprs else 'predict_{0}_with_{1}_{2}{3}.p'
    if args.biomarkers is None:
        predict_file_basename = predict_file_trunk.format(predict_biomarker_str,
                                                          method, '_'.join(visits),
                                                          '_last' if use_last_visit else '')
    else:
        estimate_biomarkers_string = '_'.join(args.biomarkers).replace(' ', '_')
        predict_file_basename = predict_file_trunk.format(predict_biomarker_str,
                                                          estimate_biomarkers_string,
                                                          '_'.join(visits),
                                                          '_last' if use_last_visit else '')
    prediction_file = os.path.join(adni.eval_folder, predict_file_basename)

    # Read if predictions exist, else recompute
    if os.path.isfile(prediction_file) and not recompute_predictions:
        # Read biomarker predictions from file
        print log.INFO, 'Reading {0} predictions from {1}...'.format(biomarker, prediction_file)
        (diagnoses, values_observed, values_naive, values_model) = pickle.load(open(prediction_file, 'rb'))
    else:
        next_visit = get_predicted_visit(visits)
        print log.INFO, 'Predicting {0} at {1}...'.format(biomarker, next_visit)

        # Get mean changes from file
        mean_changes_file = os.path.join(adni.eval_folder, 'mean_changes.p')
        mean_changes = pickle.load(open(mean_changes_file, 'rb'))

        # Get DPI estimates
        rids, diagnoses_all, dpis, dprs, _, _ = get_progress_estimates(args, visits, method=method,
                                                                       recompute_estimates=recompute_estimates,
                                                                       estimate_dprs=estimate_dprs,
                                                                       consistent_data=consistent_data)

        # Collect MMSE data for test
        data_handler = DataHandler.get_data_handler()
        model = ProgressionModel(biomarker, data_handler.get_model_file(biomarker))
        measurements = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24', 'm36'],
                                                             biomarkers=[biomarker],
                                                             select_complete=True)

        print log.INFO, 'Predicting {0} for {1}'.format(biomarker, next_visit)
        values_observed = []
        values_model = []
        values_naive = []
        diagnoses = []
        for diagnosis, rid, dpi, dpr in zip(diagnoses_all, rids, dpis, dprs):
            if rid in measurements:
                # Get real biomarker value value at next visit
                scantime_first_visit = measurements[rid][visits[0]]['scantime']
                scantime_next_visit = measurements[rid][next_visit]['scantime']
                progress_next_visit = ModelFitter.scantime_to_progress(scantime_next_visit, scantime_first_visit, dpi, dpr)
                value_observed = measurements[rid][next_visit][biomarker]
                values_observed.append(value_observed)

                # Predict biomarker value value at next visit
                if use_last_visit:
                    value = measurements[rid][visits[-1]][biomarker]
                    scantime = measurements[rid][visits[-1]]['scantime']
                    progress = ModelFitter.scantime_to_progress(scantime, scantime_first_visit, dpi, dpr)
                    mean_quantile = model.approximate_quantile(progress, value)
                else:
                    mean_quantile = 0.0
                    for visit in visits:
                        value = measurements[rid][visit][biomarker]
                        scantime = measurements[rid][visit]['scantime']
                        progress = ModelFitter.scantime_to_progress(scantime, scantime_first_visit, dpi, dpr)
                        mean_quantile += model.approximate_quantile(progress, value)
                    mean_quantile /= len(visits)

                value_model = model.get_value_at_quantile(progress_next_visit, mean_quantile)
                values_model.append(value_model)

                # Predict biomarker value naively
                if naive_use_diagnosis:
                    mean_change = mean_changes[biomarker][diagnosis]
                else:
                    mean_change = mean_changes[biomarker][0.66]

                if use_last_visit:
                    x = measurements[rid][visits[-1]]['scantime']
                    y = measurements[rid][visits[-1]][biomarker]
                    intercept = -(mean_change * x - y)
                else:
                    x = np.zeros(len(visits))
                    y = np.zeros(len(visits))
                    for i, visit in enumerate(visits):
                        x[i] = measurements[rid][visit]['scantime']
                        y[i] = measurements[rid][visit][biomarker]
                    intercept = -np.sum(mean_change * x - y) / len(x)

                value_naive = intercept + mean_change * measurements[rid][next_visit]['scantime']
                values_naive.append(value_naive)

                # Plot estimates
                if plot_predictions and diagnosis != 0.0:
                    plot_prediction_approaches(model, visits, biomarker, measurements[rid],
                                               dpi, dpr, value_model, value_naive,
                                               mean_quantile, mean_change, intercept)

                # Append diagnosis
                diagnoses.append(diagnosis)

                # Print result
                print log.RESULT, '{0} for subject {1}: Observed: {2}, Naive {3}, Model: {4}'.format(biomarker, rid, value_observed, value_naive, value_model)

        # Save results
        print log.INFO, 'Saving {0} predictions to {1}...'.format(biomarker, prediction_file)
        pickle.dump((diagnoses, values_observed, values_naive, values_model), open(prediction_file, 'wb'))

    diagnoses = np.array(diagnoses)
    values_observed = np.array(values_observed)
    values_naive = np.array(values_naive)
    values_model = np.array(values_model)

    # Exclude healthy subjects
    if exclude_cn:
        indices = np.where(diagnoses != 0.0)
        diagnoses = diagnoses[indices]
        values_observed = values_observed[indices]
        values_naive = values_naive[indices]
        values_model = values_model[indices]

    return diagnoses, values_observed, values_naive, values_model


def get_predicted_visit(visits):
    if visits[-1] == 'bl':
        return 'm12'
    elif visits[-1] == 'm12':
        return 'm24'
    elif visits[-1] == 'm24':
        return 'm36'
    else:
        print log.ERROR, 'Invalid last visit!'
        return None


def plot_prediction_approaches(model, visits, predict_biomarker, rid_measurements, dpi, dpr, value_model, value_naive, mean_quantile, change, intercept):
    next_visit = get_predicted_visit(visits)
    scantime_first_visit = rid_measurements[visits[0]]['scantime']
    scantime_next_visit = rid_measurements[next_visit]['scantime']
    progress_first_visit = ModelFitter.scantime_to_progress(scantime_first_visit, scantime_first_visit, dpi, dpr)
    progress_next_visit = ModelFitter.scantime_to_progress(scantime_next_visit, scantime_first_visit, dpi, dpr)
    progress_linspace = np.linspace(progress_first_visit - 200, progress_next_visit + 200, 100)

    fig, ax = plt.subplots()
    setup_axes(plt, ax)

    color_mapper = cm.ScalarMappable(cmap=plt.get_cmap(aplt.progression_cmap),
                                     norm=colors.Normalize(vmin=0.0, vmax=1.0))

    # Plot the percentile curves of the fitted model
    quantiles = [0.1, 0.25, 0.5, 0.75, 0.9]
    grey_values = ['0.8', '0.6', '0.4', '0.62', '0.84']
    for grey_value, quantile in zip(grey_values, quantiles):
        curve = model.get_quantile_curve(progress_linspace, quantile)
        ax.plot(progress_linspace, curve, zorder=1, color=grey_value)

    # Collect points
    progr_points = []
    value_points = []
    diagn_points = []
    for visit in visits + [next_visit]:
        value_points.append(rid_measurements[visit][predict_biomarker])
        progr_points.append(ModelFitter.scantime_to_progress(rid_measurements[visit]['scantime'], scantime_first_visit, dpi, dpr))
        diagn_points.append(rid_measurements[visit]['DX.scan'])

    # Collect lines
    predict_diagnosis = rid_measurements[next_visit]['DX.scan']
    predict_linspace = np.linspace(progress_first_visit, progress_next_visit, 50)
    curve = [model.get_value_at_quantile(p, mean_quantile) for p in predict_linspace]
    line = [change * ModelFitter.progress_to_scantime(p, scantime_first_visit, dpi, dpr) + intercept for p in predict_linspace]

    # Plot model and linear prediction line
    ax.plot(predict_linspace, line, zorder=1, linestyle='--', linewidth=2, color='k',
            label='naive prediction')
    ax.plot(predict_linspace, curve, zorder=1, linestyle='-', linewidth=2, color='k',
            label='model-based prediction')
    ax.scatter(progr_points, value_points, zorder=2, s=50.0,
               c=[color_mapper.to_rgba(d) for d in diagn_points], edgecolor='none')

    # Plot the predicted values
    ax.scatter([progress_next_visit], [value_naive], zorder=2, s=50.0, c='w',
               edgecolor=color_mapper.to_rgba(predict_diagnosis))
    ax.scatter([progress_next_visit], [value_model], zorder=2, s=50.0, c='w',
               edgecolor=color_mapper.to_rgba(predict_diagnosis))

    plt.tight_layout()
    plt.legend()
    plt.show()
    plt.close(fig)


# ------------------------------------------------------------------------------
#  Synthetic data
# ------------------------------------------------------------------------------
def generate_synth_model(args, data_handler, biomarker, num_samples=1000,
                         sampling='longitudinal', rate_sigma=0.0, conversion_sigma=0.0, run=0):
    # Get model and samples file
    model_file_experiment = data_handler.get_model_file(biomarker, num_samples=num_samples,
                                                        sampling=sampling, rate_sigma=rate_sigma,
                                                        conversion_sigma=conversion_sigma, run=run)
    samples_file_experiment = data_handler.get_samples_file(biomarker, num_samples=num_samples,
                                                            sampling=sampling, rate_sigma=rate_sigma,
                                                            conversion_sigma=conversion_sigma, run=run)

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
            conv_sigma_str = ' --conversion_sigma {0}'.format(
                conversion_sigma) if conversion_sigma is not None and conversion_sigma > 0.0 else ''

            call('{0}/vgam_generate_synth_data.py {1}{2}{3}{4}{5}'.format(
                exec_folder, biomarker_str, num_samples_str, sampling_str,
                rate_sigma_str, conv_sigma_str), shell=True)
            call('{0}/vgam_estimate_curves.py -m synth {1}'.format(
                exec_folder, biomarker_str), shell=True)

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


def generate_synth_test_data(args, biomarkers, num_test_samples, number_of_visits, run):
    biomarkers_str = '_'.join(biomarkers)
    test_data_filename = 'test_data_{0}_{1}_{2}_{3}.p'.format(biomarkers_str, num_test_samples, number_of_visits, run)
    test_data_folder = adni.make_dir(adni.eval_folder, biomarkers_str)
    test_data_file = os.path.join(test_data_folder, test_data_filename)

    if os.path.isfile(test_data_file) and not args.recompute_test_data:
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


# ------------------------------------------------------------------------------
#    Plotting
# ------------------------------------------------------------------------------
def setup_axes(plt, ax, xgrid=True, ygrid=True):
    assert isinstance(ax, mpl.axes.Axes)

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_color('none')
    if xgrid:
        ax.xaxis.grid(True, linestyle=':', which='major', color='grey', alpha=0.5)
    if ygrid:
        ax.yaxis.grid(True, linestyle=':', which='major', color='grey', alpha=0.5)
    ax.tick_params(axis='both', which='both', bottom='on', top='off', left='off', right='off')

    plt.setp(ax.get_yticklabels(), rotation=45, fontsize=11)
    plt.setp(ax.get_xticklabels(), fontsize=11)


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
