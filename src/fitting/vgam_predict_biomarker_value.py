#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
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
    parser.add_argument('--predict_biomarker', type=str, default='MMSE', help='the biomarker to predict')
    parser.add_argument('--estimate_dprs', action='store_true', help='estimate dpis and dprs')
    parser.add_argument('--recompute_dpis', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--exclude_cn', action='store_true', help='exclude healthy subjects from analysis')
    parser.add_argument('--plot_steps', type=int, default=15, help='number of steps for the DPI scale')
    parser.add_argument('--output_file', type=str, default=None, help='filename of the output file')
    args = parser.parse_args()

    diagnoses, values_observed, values_naive, values_predicted = get_biomarker_predictions(args)
    # plot_biomarker_predictions(args, diagnoses, values_observed, values_predicted)
    analyse_biomarker_predictions(args, diagnoses, values_observed, values_naive, values_predicted)


def get_biomarker_predictions(args):
    biomarker = args.predict_biomarker
    biomarker_str = biomarker.replace(' ', '_')
    predict_file_basename = 'predict_{0}_with_{1}_{2}.p'.format(biomarker_str, args.method, '_'.join(args.visits))
    predict_file = os.path.join(adni.eval_folder, predict_file_basename)

    mean_changes_file = os.path.join(adni.eval_folder, 'mean_changes.p')
    mean_changes = pickle.load(open(mean_changes_file, 'rb'))

    if os.path.isfile(predict_file) and not args.recompute_dpis:
        # Read test results from file
        print log.INFO, 'Reading {0} predictions from {1}...'.format(biomarker, predict_file)
        (diagnoses, values_observed, values_naive, values_predicted) = pickle.load(open(predict_file, 'rb'))
    else:
        # Get DPI estimates
        rids, diagnoses_all, dpis, _, _ = ve.get_dpi_estimates(args)

        # Collect MMSE data for test
        data_handler = DataHandler.get_data_handler()
        next_visit = get_next_visit(args.visits)
        predict_model = ProgressionModel(args.predict_biomarker,
                                         data_handler.get_model_file(biomarker))
        predict_measurements = data_handler.get_measurements_as_dict(visits=args.visits + [next_visit],
                                                                     biomarkers=[biomarker],
                                                                     select_complete=True)

        print log.INFO, 'Predicting {0} for {1}'.format(biomarker, next_visit)
        values_observed = []
        values_predicted = []
        values_naive = []
        diagnoses = []
        for diagnosis, rid, dpi in zip(diagnoses_all, rids, dpis):
            if rid in predict_measurements:
                # Get real biomarker value value at next visit
                progression_observed = dpi + predict_measurements[rid][next_visit]['scantime']
                value_observed = predict_measurements[rid][next_visit][biomarker]
                values_observed.append(value_observed)

                # Predict biomarker value value at next visit
                mean_quantile = 0.0
                for visit in args.visits:
                    value = predict_measurements[rid][visit][biomarker]
                    progression = dpi + predict_measurements[rid][visit]['scantime']
                    mean_quantile += predict_model.approximate_quantile(progression, value)
                mean_quantile /= len(args.visits)

                value_predicted = predict_model.get_value_at_quantile(progression_observed, mean_quantile)
                values_predicted.append(value_predicted)

                # Predict biomarker value naively
                num_visits = len(args.visits)
                x = np.zeros(num_visits)
                y = np.zeros(num_visits)
                for i, visit in enumerate(args.visits):
                    x[i] = predict_measurements[rid][visit]['scantime']
                    y[i] = predict_measurements[rid][visit][biomarker]
                mean_change = mean_changes[biomarker][diagnosis]
                y = y + mean_change * x
                A = np.array([np.zeros(num_visits), np.ones(num_visits)])
                intercept = np.linalg.lstsq(A.T, y)[0][1]

                value_naive = intercept + mean_change * progression_observed
                values_naive.append(value_naive)

                # Append diagnosis
                diagnoses.append(diagnosis)

                # Print result
                print log.RESULT, '{0} for subject {1}: Observed: {2}, Naive {3}, Predicted: {4}'.format(biomarker, rid, value_observed, value_naive, value_predicted)

        # Save results
        print log.INFO, 'Saving {0} predictions to {1}...'.format(biomarker, predict_file)
        pickle.dump((diagnoses, values_observed, values_naive, values_predicted), open(predict_file, 'wb'))

    return diagnoses, values_observed, values_naive, values_predicted


def plot_biomarker_predictions(args, diagnoses, values_observed, values_predicted):
    # Setup plot
    fig, ax = plt.subplots()
    ve.setup_axes(plt, ax)

    ax.set_title('{0} estimation based on {1} with {2} visits'.format(args.predict_biomarker,
                                                                      args.method,
                                                                      len(args.visits)))
    ax.set_xlabel('Observed {0}'.format(args.predict_biomarker))
    ax.set_ylabel('Predicted {0}'.format(args.predict_biomarker))

    # Draw plot
    plt.scatter(values_observed, values_predicted, s=15.0,
                c=diagnoses, edgecolor='none', vmin=0.0, vmax=1.0,
                cmap=aplt.progression_cmap, alpha=0.5)

    # Show or save the plot
    plt.tight_layout()
    if args.output_file is not None:
        plt.savefig(args.output_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


def analyse_biomarker_predictions(args, diagnoses, values_observed, values_naive, values_predicted):
    values_observed = np.array(values_observed)
    values_naive = np.array(values_naive)
    values_predicted = np.array(values_predicted)

    # Exclude health subjects
    if args.exclude_cn:
        indices = np.where(np.array(diagnoses) != 0.0)
        values_observed = values_observed[indices]
        values_naive = values_naive[indices]
        values_predicted = values_predicted[indices]

    print log.RESULT, '-----------------------'
    print log.RESULT, 'Naive prediction:'
    mean_error = np.sum(np.abs(values_observed - values_naive)) / len(values_observed)
    rms_error = math.sqrt(np.sum(np.square(values_observed - values_naive)) / len(values_observed))
    print log.RESULT, 'Mean error: {0}'.format(mean_error)
    print log.RESULT, 'RMS error: {0}'.format(rms_error)

    personr = stats.stats.pearsonr(values_observed, values_naive)
    print log.RESULT, 'Pearsons correlation coefficient: {0}'.format(personr[0])
    print log.RESULT, '2-tailed p-value: {0}'.format(personr[1])

    print log.RESULT, '-----------------------'
    print log.RESULT, 'Model-based prediction:'
    mean_error = np.sum(np.abs(values_observed - values_predicted)) / len(values_observed)
    rms_error = math.sqrt(np.sum(np.square(values_observed - values_predicted)) / len(values_observed))
    print log.RESULT, 'Mean error: {0}'.format(mean_error)
    print log.RESULT, 'RMS error: {0}'.format(rms_error)

    personr = stats.stats.pearsonr(values_observed, values_predicted)
    print log.RESULT, 'Pearsons correlation coefficient: {0}'.format(personr[0])
    print log.RESULT, '2-tailed p-value: {0}'.format(personr[1])


def get_next_visit(visits):
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
