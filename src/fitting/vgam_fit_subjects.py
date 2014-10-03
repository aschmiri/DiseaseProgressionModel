#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
from common import log as log
from common import adni_plot as aplt
from common import adni_tools as adni
from vgam.datahandler import DataHandler
from vgam.progressionmodel import MultiBiomarkerProgressionModel
from vgam.modelfitter import ModelFitter


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('visits', nargs='+', type=str, help='the viscodes to be sampled')
    parser.add_argument('--recompute_dpis', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--plot_steps', type=int, default=15, help='number of steps for the DPI scale')
    parser.add_argument('--output_file', type=str, default=None, help='filename of the output file')
    args = parser.parse_args()

    evaluation_basename = 'eval_{0}__{1}.p'.format(args.method, '_'.join(args.visits))
    evaluation_file = os.path.join(adni.eval_folder, evaluation_basename)

    if os.path.isfile(evaluation_file) and not args.recompute_dpis:
        # Read test results from file
        print log.INFO, 'Reading test results from {0}...'.format(evaluation_file)
        (dpis, diagnoses, mean_min, mean_max) = pickle.load(open(evaluation_file, 'rb'))
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

        # Calculate data
        mean_min = model.get_mean_min_progression()
        mean_max = model.get_mean_max_progression()
        dpis, diagnoses = main_estimate_dpi(measurements, args.visits, fitter)

        # Save results
        print log.INFO, 'Saving test results to {0}...'.format(evaluation_file)
        pickle.dump((dpis, diagnoses, mean_min, mean_max), open(evaluation_file, 'wb'))

    analyse_dpis(args, dpis, diagnoses, mean_min, mean_max)

    # Get RCDs
    # rcds = data_handler.get_mmse_decline_as_dict()
    # main_estimate_dpi_dpr(measurements, viscodes, fitter, rcds)
    # main_evaluate_scalings(measurements, viscodes, fitter, rcds)


def main_estimate_dpi(measurements, viscodes, fitter):
    # Test all available subjects
    dpis = []
    diagnoses = []
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

        print log.RESULT, 'Estimated DPI: {0}, Diagnosis: {1}'.format(dpi, diagnosis)
        if dpi is not None:
            dpis.append(dpi)
            diagnoses.append(diagnosis)

    return dpis, diagnoses


def main_estimate_dpi_dpr(measurements, viscodes, fitter, rcds):
    # Test all available subjects
    dpis = []
    dprs = []
    progresses = []
    rcdnum = []
    for rid in measurements:
        print log.INFO, 'Estimating DPI and DPR for subject {0}...'.format(rid)
        try:
            progress = measurements[rid]['bl']['progress']
            rcd = rcds[rid]
        except KeyError:
            continue

        if not set(viscodes).issubset(set(measurements[rid].keys())):
            print log.WARNING, 'Not all viscodes {0} available for subject {1}!'.format(viscodes, rid)
            continue

        samples = {}
        for viscode in viscodes:
            samples.update({viscode: measurements[rid][viscode]})
        dpi, dpr = fitter.get_dpi_dpr_for_samples(samples)

        print log.RESULT, 'Estimated DPI: {0}, DPR: {1}, Progress: {2}'.format(dpi, dpr, progress)
        if dpi is not None and dpr is not None:
            dpis.append(dpi)
            dprs.append(dpr)
            progresses.append(progress)
            rcdnum.append(rcd)

    rms_error = np.sqrt(np.sum(np.square(np.array(progresses) - np.array(dpis))) / len(dpis))
    mean_error = np.sum(np.abs(np.array(progresses) - np.array(dpis))) / len(dpis)
    print log.RESULT, 'Mean error: {0}, RMS error: {1}'.format(mean_error, rms_error)

    # Plot the results
    plot_rcds(dpis, dprs, rcdnum)


def main_evaluate_scalings(measurements, fitter, rcds):
    measurements = fitter.get_scaled_measurements(measurements)

    # Test all available subjects
    scalings = []
    rcdnum = []
    for rid in measurements:
        if rid in rcds:
            scalings.append(measurements[rid]['bl']['scaling'])
            rcdnum.append(rcds[rid])

    # Plot the results
    plt.title('Correlation between scaling value and MMS decline')
    plt.xlabel('Estimated scaling values')
    plt.ylabel('MMSE decline')
    plt.scatter(scalings, rcdnum, s=50, linewidths=0, alpha=0.5)
    plt.show()


def analyse_dpis(args, dpis, diagnoses, mean_min, mean_max):
    dpi_range = float(ModelFitter.TEST_DPI_MAX - ModelFitter.TEST_DPI_MIN)

    # Setup plot
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.set_title('DPI estimation using {0} at {1}'.format(args.method, ', '.join(args.visits)))
    ax.set_xlabel('DPI')
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(['CN', 'EMCI', 'LMCI', 'AD'])

    xticks = np.linspace(0, args.plot_steps, 7)
    ax.set_xticks(xticks)
    ax.set_xticklabels([int(float(tick) / args.plot_steps * dpi_range + ModelFitter.TEST_DPI_MIN) for tick in xticks])

    # Compute matrix
    diagnosis_indices = {0.0: 0, 0.25: 1, 0.5: 1, 0.75: 2, 1.0: 3}
    matrix = np.zeros((4, args.plot_steps + 1))
    mean = np.zeros(4)
    for dpi, diag in zip(dpis, diagnoses):
        diagnosis_index = diagnosis_indices[diag]
        dpi_index = round((dpi - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps)
        matrix[diagnosis_index, dpi_index] += 1.0
        mean[diagnosis_index] += dpi

    # Draw annotations
    dpis = np.array(dpis)
    diagnoses = np.array(diagnoses)
    means = []
    stds = []
    for diag in [0.0, 0.25, 0.75, 1.0]:
        i = diagnosis_indices[diag]
        num_subjects = np.sum(matrix[i])
        matrix[i] /= num_subjects

        indices = np.where(diagnoses == diag)
        mean = (np.mean(dpis[indices]) - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps
        std = np.std(dpis[indices]) / dpi_range * args.plot_steps
        means.append(mean)
        stds.append(std)
        # ax.text(args.plot_steps + 1, i, '$n={0}$'.format(int(num_subjects)))

    plt.errorbar(means, [0, 1, 2, 3], xerr=stds, fmt=None, ecolor='r', elinewidth=2, capsize=4, capthick=2)
    plt.plot(means, [0, 1, 2, 3], linestyle='', color='r', marker='|', markersize=15, markeredgewidth=2)
    plt.axvline((mean_min - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps, color='k', linestyle=':', alpha=0.6)
    plt.axvline((mean_max - ModelFitter.TEST_DPI_MIN) / dpi_range * args.plot_steps, color='k', linestyle=':', alpha=0.6)
    plt.axvline(-ModelFitter.TEST_DPI_MIN / dpi_range * args.plot_steps, color='k', linestyle='-', alpha=0.6)
    plt.imshow(matrix, cmap=plt.get_cmap('Greys'), interpolation='nearest')

    # Draw or save the plot
    plt.tight_layout()
    if args.output_file is not None:
        plt.savefig(args.output_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)

    # Calculate accuracy
    misclassified = 0
    for dpi, diagnosis in zip(dpis, diagnoses):
        if dpi > 0 and diagnosis == 0.5 or dpi <= 0 and diagnosis == 1.0:
            misclassified += 1
    print log.RESULT, 'Accuracy: {0}'.format(1.0 - float(misclassified) / float(len(dpis)))


def plot_rcds(dpis, dprs, rcds):
    plt.title('DPI and DPR with RCD')
    plt.xlabel('Estimated DPI')
    plt.ylabel('Estimated DPR')
    plt.scatter(dpis, dprs, c=rcds, cmap=aplt.progression_cmap, s=50, linewidths=0)
    plt.show()


if __name__ == '__main__':
    main()
