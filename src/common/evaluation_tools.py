#! /usr/bin/env python2.7
import os.path
import pickle
import numpy as np
from common import log as log
from common.datahandler import DataHandler
from common.progressionmodel import ProgressionModel
from common.progressionmodel import MultiBiomarkerProgressionModel
from common.modelfitter import ModelFitter


# ------------------------------------------------------------------------------
#  Progression estimation
# ------------------------------------------------------------------------------
def get_progress_estimates(visits,
                           method=None, biomarkers=None, phase=None,
                           recompute_estimates=False,
                           estimate_dprs=False, consistent_data=False,
                           select_training_set=False, select_test_set=False):
    # Get data handler and biomarker names
    data_handler = DataHandler.get_data_handler(method=method,
                                                biomarkers=biomarkers,
                                                phase=phase)

    # Get filename
    estimates_file_trunk = 'estimate_dpi_dpr_with_{0}_{1}.p' if estimate_dprs else 'estimate_dpi_with_{0}_{1}.p'
    if biomarkers is None:
        estimates_file_basename = estimates_file_trunk.format(method, '_'.join(visits))
    else:
        biomarkers_string = '_'.join(biomarkers).replace(' ', '_')
        estimates_file_basename = estimates_file_trunk.format(biomarkers_string, '_'.join(visits))
    estimates_file = os.path.join(data_handler.get_eval_folder(), estimates_file_basename)

    # Read if estimates exist, else recompute
    if os.path.isfile(estimates_file) and not recompute_estimates:
        # Read test results from file
        print log.INFO, 'Reading DPI{0} estimations from {1}...'.format('\DPR' if estimate_dprs else '', estimates_file)
        (rids, diagnoses, dpis, dprs, mean_min, mean_max) = pickle.load(open(estimates_file, 'rb'))
    else:
        # Collect data for test
        biomarkers = data_handler.get_biomarker_names()
        measurements = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24'],
                                                             biomarkers=biomarkers,
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
            rids, diagnoses, dpis = estimate_dpis(measurements, visits, fitter, phase=phase)
            dprs = np.ones(len(dpis)).tolist()
        else:
            rids, diagnoses, dpis, dprs = estimate_dpis_dprs(measurements, visits, fitter, phase=phase)

        print log.INFO, 'Saving DPI{0} estimations to {1}...'.format('\DPR' if estimate_dprs else '', estimates_file)
        pickle.dump((rids, diagnoses, dpis, dprs, mean_min, mean_max), open(estimates_file, 'wb'))

    # Reduce to consistent data sets with bl, m12 and m24 samples
    if consistent_data or select_training_set or select_test_set:
        consistent_method = 'all' if consistent_data else method
        consistent_data_handler = DataHandler.get_data_handler(method=consistent_method)
        consistent_measurements = consistent_data_handler.get_measurements_as_dict(
            visits=['bl', 'm12', 'm24'],
            select_training_set=select_training_set,
            select_test_set=select_test_set,
            select_complete=True,
            no_regression=True)

        consistent_rids = []
        consistent_diagnoses = []
        consistent_dpis = []
        consistent_dprs = []
        for i, rid in enumerate(rids):
            if rid in consistent_measurements:
                consistent_rids.append(rid)
                consistent_diagnoses.append(diagnoses[i])
                consistent_dpis.append(dpis[i])
                consistent_dprs.append(dprs[i])
        rids = consistent_rids
        diagnoses = consistent_diagnoses
        dpis = consistent_dpis
        dprs = consistent_dprs

        print log.RESULT, 'Selected {0} consistent subjects.'.format(len(dpis))

    # Return results
    return rids, diagnoses, dpis, dprs, mean_min, mean_max


def estimate_dpis(measurements, viscodes, fitter, phase=None):
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
        dpi = fitter.get_dpi_for_samples(samples, phase=phase)

        print log.RESULT, 'Subject {0}: Estimated DPI: {1}, Diagnosis: {2}'.format(rid, dpi, diagnosis)
        if dpi is not None:
            rids.append(rid)
            diagnoses.append(diagnosis)
            dpis.append(dpi)

    return rids, diagnoses, dpis


def estimate_dpis_dprs(measurements, viscodes, fitter, phase=None):
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
        dpi, dpr = fitter.get_dpi_dpr_for_samples(samples, phase=phase)

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
def get_biomarker_predictions(visits, predict_biomarker,
                              method=None, biomarkers=None, phase=None,
                              recompute_estimates=False, recompute_predictions=False,
                              estimate_dprs=False, consistent_data=False, exclude_cn=False,
                              use_last_visit=False, naive_use_diagnosis=False):

    # Get prediction file
    data_handler = DataHandler.get_data_handler(method=method,
                                                biomarkers=biomarkers,
                                                phase=phase)
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
    prediction_file = os.path.join(data_handler.get_eval_folder(), predict_file_basename)

    # Read if predictions exist, else recompute
    if os.path.isfile(prediction_file) and not recompute_predictions:
        # Read biomarker predictions from file
        print log.INFO, 'Reading {0} predictions from {1}...'.format(predict_biomarker, prediction_file)
        (rids, diagnoses, values_observed, values_naive, values_model) = pickle.load(open(prediction_file, 'rb'))
    else:
        predict_visit = get_predicted_visit(visits)
        print log.INFO, 'Predicting {0} at {1}...'.format(predict_biomarker, predict_visit)

        # Get mean changes from file
        mean_changes_file = os.path.join(data_handler.get_eval_folder(), 'mean_changes.p')
        if not os.path.isfile(mean_changes_file):
            print log.ERROR, 'Mean changes unknown, run misc/compute_mean_biomarker_changes.py first!'
        mean_changes = pickle.load(open(mean_changes_file, 'rb'))

        # Get DPI estimates
        rids_all, diagnoses_all, dpis, dprs, _, _ = get_progress_estimates(visits,
                                                                           method=method,
                                                                           biomarkers=biomarkers,
                                                                           phase=phase,
                                                                           recompute_estimates=recompute_estimates,
                                                                           estimate_dprs=estimate_dprs,
                                                                           consistent_data=consistent_data)

        # Collect biomarker data for test
        measurements = data_handler.get_measurements_as_dict(visits=visits + [predict_visit],
                                                             biomarkers=[predict_biomarker],
                                                             select_complete=True)
        model = ProgressionModel(predict_biomarker, data_handler.get_model_file(predict_biomarker))

        print log.INFO, 'Predicting {0} for {1}'.format(predict_biomarker, predict_visit)
        rids = []
        diagnoses = []
        values_observed = []
        values_model = []
        values_naive = []
        for rid, diagnosis, dpi, dpr in zip(rids_all, diagnoses_all, dpis, dprs):
            if rid in measurements:
                # Get real biomarker value value at next visit
                scantime_first_visit = measurements[rid][visits[0]]['scantime']
                scantime_next_visit = measurements[rid][predict_visit]['scantime']
                progress_next_visit = ModelFitter.scantime_to_progress(scantime_next_visit, scantime_first_visit, dpi, dpr)
                value_observed = measurements[rid][predict_visit][predict_biomarker]
                values_observed.append(value_observed)

                # Predict biomarker value value at next visit
                if use_last_visit:
                    value = measurements[rid][visits[-1]][predict_biomarker]
                    scantime = measurements[rid][visits[-1]]['scantime']
                    progress = ModelFitter.scantime_to_progress(scantime, scantime_first_visit, dpi, dpr)
                    mean_quantile = model.approximate_quantile(progress, value)
                else:
                    mean_quantile = 0.0
                    for visit in visits:
                        value = measurements[rid][visit][predict_biomarker]
                        scantime = measurements[rid][visit]['scantime']
                        progress = ModelFitter.scantime_to_progress(scantime, scantime_first_visit, dpi, dpr)
                        mean_quantile += model.approximate_quantile(progress, value)
                    mean_quantile /= len(visits)

                value_model = model.get_value_at_quantile(progress_next_visit, mean_quantile)
                values_model.append(value_model)

                # Predict biomarker value naively
                if naive_use_diagnosis:
                    mean_change = mean_changes[predict_biomarker][diagnosis]
                else:
                    mean_change = mean_changes[predict_biomarker][0.66]

                if use_last_visit:
                    x = measurements[rid][visits[-1]]['scantime']
                    y = measurements[rid][visits[-1]][predict_biomarker]
                    intercept = -(mean_change * x - y)
                else:
                    x = np.zeros(len(visits))
                    y = np.zeros(len(visits))
                    for i, visit in enumerate(visits):
                        x[i] = measurements[rid][visit]['scantime']
                        y[i] = measurements[rid][visit][predict_biomarker]
                    intercept = -np.sum(mean_change * x - y) / len(x)

                value_naive = intercept + mean_change * measurements[rid][predict_visit]['scantime']
                values_naive.append(value_naive)

                # Append rid and diagnosis
                rids.append(rid)
                diagnoses.append(diagnosis)

                # Print result
                print log.RESULT, '{0} for subject {1}: Observed: {2}, Naive {3}, Model: {4}'.format(predict_biomarker, rid, value_observed, value_naive, value_model)

        # Save results
        print log.INFO, 'Saving {0} predictions to {1}...'.format(predict_biomarker, prediction_file)
        pickle.dump((rids, diagnoses, values_observed, values_naive, values_model), open(prediction_file, 'wb'))

    rids = np.array(rids)
    diagnoses = np.array(diagnoses)
    values_observed = np.array(values_observed)
    values_naive = np.array(values_naive)
    values_model = np.array(values_model)

    # Exclude healthy subjects
    if exclude_cn:
        indices = np.where(diagnoses != 0.0)
        rids = rids[indices]
        diagnoses = diagnoses[indices]
        values_observed = values_observed[indices]
        values_naive = values_naive[indices]
        values_model = values_model[indices]

    return rids, diagnoses, values_observed, values_naive, values_model


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
