#! /usr/bin/env python2.7
import os.path
import argparse
import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from common import evaluation_tools as et
from common import plotting_tools as pt
from common import log as log
from common.datahandler import DataHandler


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('visits', nargs='+', type=str, help='the viscodes of the visits that are available')
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-p', '--phase', default=None, choices=DataHandler.get_phase_choices(), help='the phase for which the model is to be trained')
    parser.add_argument('--predict_biomarker', type=str, default='MMSE', help='the biomarker to predict')
    parser.add_argument('--recompute_estimates', action='store_true', help='recompute the dpi / dpr estimations')
    parser.add_argument('--recompute_predictions', action='store_true', help='recompute the biomarker predictions')
    parser.add_argument('--estimate_dprs', action='store_true', help='estimate dpis and dprs')
    parser.add_argument('--consistent_data', action='store_true', help='use only subjects with bl, m12 and m24 visits')
    parser.add_argument('--exclude_cn', action='store_true', help='exclude healthy subjects from analysis')
    parser.add_argument('--use_last_visit', action='store_true', help='use only the last visit for prediction')
    parser.add_argument('--naive_use_diagnosis', action='store_true', help='use the specific mean change for the diagnosis')
    parser.add_argument('--no_plot', action='store_true', help='do not plot the results')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output file')
    parser.add_argument('--latex_file', type=str, default=None, help='add output to a LaTeX file')
    args = parser.parse_args()

    _, diagnoses, values_observed, values_naive, values_model = \
        et.get_biomarker_predictions(args.visits, args.predict_biomarker,
                                     method=args.method,
                                     biomarkers=args.biomarkers,
                                     phase=args.phase,
                                     recompute_estimates=args.recompute_estimates,
                                     recompute_predictions=args.recompute_predictions,
                                     estimate_dprs=args.estimate_dprs,
                                     select_test_set=True,
                                     consistent_data=args.consistent_data,
                                     exclude_cn=args.exclude_cn,
                                     use_last_visit=args.use_last_visit,
                                     naive_use_diagnosis=args.naive_use_diagnosis)
    if not args.no_plot:
        plot_biomarker_predictions(args, diagnoses, values_observed, values_model)
    analyse_biomarker_predictions(args, diagnoses, values_observed, values_naive, values_model)


def plot_biomarker_predictions(args, diagnoses, values_observed, values_model):
    # Setup plot
    fig, ax = plt.subplots()
    pt.setup_axes(plt, ax)

    ax.set_title('{0} estimation based on {1} with {2} visits'.format(args.predict_biomarker,
                                                                      args.method,
                                                                      len(args.visits)))
    ax.set_xlabel('Observed {0}'.format(args.predict_biomarker))
    ax.set_ylabel('Predicted {0}'.format(args.predict_biomarker))

    # Draw plot
    plt.scatter(values_observed, values_model, s=15.0,
                c=diagnoses, edgecolor='none', vmin=0.0, vmax=1.0,
                cmap=pt.progression_cmap, alpha=0.5)

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
    data_handler = DataHandler.get_data_handler(method=args.method,
                                                biomarkers=args.biomarkers,
                                                phase=args.phase)
    filename = os.path.join(data_handler.get_eval_folder(), args.latex_file)
    with open(filename, 'a') as latex_file:
        latex_file.write('{0} & {1} {2} & ${3:.2f}\pm{4:.2f}$ & ${5:.2f}$ & ${6:.2f}\pm{7:.2f}$ & ${8:.2f}$ & {9}\\\\\n'.format(
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


if __name__ == '__main__':
    main()
