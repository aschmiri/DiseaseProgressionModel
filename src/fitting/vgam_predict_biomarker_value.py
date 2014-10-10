#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from common import log as log
from common import adni_tools as adni
from common import adni_plot as aplt
from vgam.datahandler import DataHandler
from vgam.progressionmodel import ProgressionModel
import fitting.vgam_evaluation as ve


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('visits', nargs='+', type=str, help='the viscodes of the visits that are available')
    parser.add_argument('--predict_biomarker', type=str, default='MMSE', choices=adni.biomarker_names, help='the biomarker to predict')
    parser.add_argument('--estimate_dprs', action='store_true', help='estimate dpis and dprs')
    parser.add_argument('--recompute_estimates', action='store_true', help='recompute the dpi / dpr estimations')
    parser.add_argument('--recompute_predictions', action='store_true', help='recompute the biomarker predictions')
    parser.add_argument('--exclude_cn', action='store_true', help='exclude healthy subjects from analysis')
    parser.add_argument('--naive_use_diagnosis', action='store_true', help='use the specific mean change for the diagnosis')
    parser.add_argument('--plot_predictions', action='store_true', help='visualise the predictions of each subject')
    parser.add_argument('--no_plot', action='store_true', help='do not plot the results')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output file')
    parser.add_argument('--latex_file', type=str, default=None, help='add output to a LaTeX file')
    args = parser.parse_args()

    diagnoses, values_observed, values_naive, values_model = get_biomarker_predictions(args)
    if not args.no_plot:
        plot_biomarker_predictions(args, diagnoses, values_observed, values_model)
    analyse_biomarker_predictions(args, diagnoses, values_observed, values_naive, values_model)


def get_biomarker_predictions(args):
    biomarker = args.predict_biomarker

    # Get prediction file
    predict_biomarker_str = biomarker.replace(' ', '_')
    predict_file_trunk = 'predict_{0}_with_dpr_{1}_{2}.p' if args.estimate_dprs else 'predict_{0}_with_{1}_{2}.p'
    if args.biomarkers_name is None:
        predict_file_basename = predict_file_trunk.format(predict_biomarker_str, args.method, '_'.join(args.visits))
    else:
        estimate_biomarkers_string = '_'.join(args.biomarkers_name).replace(' ', '_')
        predict_file_basename = predict_file_trunk.format(predict_biomarker_str, estimate_biomarkers_string, '_'.join(args.visits))
    predict_file = os.path.join(adni.eval_folder, predict_file_basename)

    # Read if predictions exist, else recompute
    if os.path.isfile(predict_file) and not args.recompute_predictions:
        # Read biomarker predictions from file
        print log.INFO, 'Reading {0} predictions from {1}...'.format(biomarker, predict_file)
        (diagnoses, values_observed, values_naive, values_model) = pickle.load(open(predict_file, 'rb'))
    else:
        next_visit = get_predicted_visit(args.visits)
        print log.INFO, 'Predicting {0} at {1}...'.format(biomarker, next_visit)

        # Get mean changes from file
        mean_changes_file = os.path.join(adni.eval_folder, 'mean_changes.p')
        mean_changes = pickle.load(open(mean_changes_file, 'rb'))

        # Get DPI estimates
        rids, diagnoses_all, dpis, _, _, _ = ve.get_dpi_estimates(args)

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
        for diagnosis, rid, dpi in zip(diagnoses_all, rids, dpis):
            if rid in measurements:
                # Get real biomarker value value at next visit
                progression_observed = dpi + measurements[rid][next_visit]['scantime']
                value_observed = measurements[rid][next_visit][biomarker]
                values_observed.append(value_observed)

                # Predict biomarker value value at next visit
                mean_quantile = 0.0
                for visit in args.visits:
                    value = measurements[rid][visit][biomarker]
                    progression = dpi + measurements[rid][visit]['scantime']
                    mean_quantile += model.approximate_quantile(progression, value)
                mean_quantile /= len(args.visits)

                value_model = model.get_value_at_quantile(progression_observed, mean_quantile)
                values_model.append(value_model)

                # Predict biomarker value naively
                num_visits = len(args.visits)
                x = np.zeros(num_visits)
                y = np.zeros(num_visits)
                for i, visit in enumerate(args.visits):
                    x[i] = measurements[rid][visit]['scantime']
                    y[i] = measurements[rid][visit][biomarker]
                if args.naive_use_diagnosis:
                    mean_change = mean_changes[biomarker][diagnosis]
                else:
                    mean_change = mean_changes[biomarker][0.66]
                intercept = -np.sum(mean_change * x - y) / len(x)

                value_naive = intercept + mean_change * measurements[rid][next_visit]['scantime']
                values_naive.append(value_naive)

                # Plot estimates
                if args.plot_predictions and diagnosis != 0.0:
                    plot_predictions(args, model, measurements[rid], dpi,
                                     value_model, value_naive,
                                     mean_quantile, mean_change, intercept)

                # Append diagnosis
                diagnoses.append(diagnosis)

                # Print result
                print log.RESULT, '{0} for subject {1}: Observed: {2}, Naive {3}, Model: {4}'.format(biomarker, rid, value_observed, value_naive, value_model)

        # Save results
        print log.INFO, 'Saving {0} predictions to {1}...'.format(biomarker, predict_file)
        pickle.dump((diagnoses, values_observed, values_naive, values_model), open(predict_file, 'wb'))

    diagnoses = np.array(diagnoses)
    values_observed = np.array(values_observed)
    values_naive = np.array(values_naive)
    values_model = np.array(values_model)

    # Exclude healthy subjects
    if args.exclude_cn:
        indices = np.where(diagnoses != 0.0)
        diagnoses = diagnoses[indices]
        values_observed = values_observed[indices]
        values_naive = values_naive[indices]
        values_model = values_model[indices]

    return diagnoses, values_observed, values_naive, values_model


def plot_biomarker_predictions(args, diagnoses, values_observed, values_model):
    # Setup plot
    fig, ax = plt.subplots()
    ve.setup_axes(plt, ax)

    ax.set_title('{0} estimation based on {1} with {2} visits'.format(args.predict_biomarker,
                                                                      args.method,
                                                                      len(args.visits)))
    ax.set_xlabel('Observed {0}'.format(args.predict_biomarker))
    ax.set_ylabel('Predicted {0}'.format(args.predict_biomarker))

    # Draw plot
    plt.scatter(values_observed, values_model, s=15.0,
                c=diagnoses, edgecolor='none', vmin=0.0, vmax=1.0,
                cmap=aplt.progression_cmap, alpha=0.5)

    # Show or save the plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def analyse_biomarker_predictions(args, diagnoses, values_observed, values_naive, values_model):
    # Compare observed with naive prediction
    print log.RESULT, '-----------------------'
    print log.RESULT, 'Naive prediction (n={0}):'.format(len(values_observed))
    results_naive = analyse_errors(values_observed, values_naive)

    # Compare observed with model-based prediction
    print log.RESULT, '-----------------------'
    print log.RESULT, 'Model-based prediction (n={0}):'.format(len(values_observed))
    results_model = analyse_errors(values_observed, values_model)

    # Print results
    if args.latex_file is not None:
        print_to_latex(args, results_naive, results_model, len(values_observed))


def analyse_errors(values1, values2):
    # Compute mean, standard deviation and RMS
    errors = np.abs(values1 - values2)
    mean_error = np.mean(errors)
    std_error = np.std(errors)
    rms_error = math.sqrt(np.mean(np.square(errors)))
    print log.RESULT, 'Mean error: {0}'.format(mean_error)
    print log.RESULT, 'Standard deviation: {0}'.format(std_error)
    print log.RESULT, 'RMS error: {0}'.format(rms_error)

    # Compute correlation coefficient
    personr = stats.stats.pearsonr(values1, values2)
    print log.RESULT, 'Pearsons correlation coefficient: {0}'.format(personr[0])
    print log.RESULT, '2-tailed p-value: {0}'.format(personr[1])

    # Return results
    results = {}
    results.update({'MEAN': mean_error})
    results.update({'STD': std_error})
    results.update({'RMS': rms_error})
    results.update({'CORR': personr[0]})
    results.update({'P': personr[1]})
    return results


def print_to_latex(args, results_naive, results_model, num_subjects):
    filename = os.path.join(adni.eval_folder, args.latex_file)
    with open(filename, 'a') as latex_file:
        latex_file.write('{0} & {1} {2} & ${3:.2f}\pm{4:.2f}$ & ${5:.2f}$ & ${6:.2f}\pm{7:.2f}$ & ${8:.2f}$ & {9:.2f}\\\\\n'.format(
                         args.predict_biomarker,
                         args.method,
                         len(args.visits),
                         results_naive['MEAN'],
                         results_naive['STD'],
                         results_naive['CORR'],
                         results_model['MEAN'],
                         results_model['STD'],
                         results_model['CORR'],
                         num_subjects))


def plot_predictions(args, model, rid_measurements, progression, value_model, value_naive, mean_quantile, change, intercept):
    predicted_visit = get_predicted_visit(args.visits)
    min_progression = progression + rid_measurements['bl']['scantime']
    max_progression = progression + rid_measurements[predicted_visit]['scantime']
    progression_linspace = np.linspace(min_progression - 200, max_progression + 200, 100)

    fig, ax = plt.subplots()
    ve.setup_axes(plt, ax)

    color_mapper = cm.ScalarMappable(cmap=plt.get_cmap(aplt.progression_cmap),
                                     norm=colors.Normalize(vmin=0.0, vmax=1.0))

    # Plot the percentile curves of the fitted model
    quantiles = [0.1, 0.25, 0.5, 0.75, 0.9]
    grey_values = ['0.8', '0.6', '0.4', '0.62', '0.84']
    for grey_value, quantile in zip(grey_values, quantiles):
        curve = model.get_quantile_curve(progression_linspace, quantile)
        ax.plot(progression_linspace, curve, zorder=1, color=grey_value)

    # Collect points
    progr_points = []
    value_points = []
    diagn_points = []
    for visit in args.visits + [predicted_visit]:
        value_points.append(rid_measurements[visit][args.predict_biomarker])
        progr_points.append(rid_measurements[visit]['scantime'] + progression)
        diagn_points.append(rid_measurements[visit]['DX.scan'])

    # Collect lines
    predict_diagnosis = rid_measurements[predicted_visit]['DX.scan']
    predict_progression = progression + rid_measurements[predicted_visit]['scantime']
    predict_linspace = np.linspace(progr_points[0], predict_progression, 50)
    curve = [model.get_value_at_quantile(p, mean_quantile) for p in predict_linspace]
    line = [change * (p - progression) + intercept for p in predict_linspace]

    # Plot model and linear prediction line
    ax.plot(predict_linspace, line, zorder=1, linestyle='--', linewidth=2, color='k',
            label='naive prediction')
    ax.plot(predict_linspace, curve, zorder=1, linestyle='-', linewidth=2, color='k',
            label='model-based prediction')
    ax.scatter(progr_points, value_points, zorder=2, s=50.0,
               c=[color_mapper.to_rgba(d) for d in diagn_points], edgecolor='none')

    # Plot the predicted values
    ax.scatter([predict_progression], [value_naive], zorder=2, s=50.0, c='w',
               edgecolor=color_mapper.to_rgba(predict_diagnosis))
    ax.scatter([predict_progression], [value_model], zorder=2, s=50.0, c='w',
               edgecolor=color_mapper.to_rgba(predict_diagnosis))

    plt.tight_layout()
    plt.legend()
    plt.show()
    plt.close(fig)


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


if __name__ == '__main__':
    main()
