#! /usr/bin/env python2.7
import argparse
import numpy as np
import matplotlib.pyplot as plt
from common import log as log
from common import adni_plot as aplt
from vgam.datahandler import DataHandler
from vgam.progressionmodel import MultiBiomarkerProgressionModel
from vgam.modelfitter import ModelFitter


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('visits', nargs='+', type=str, help='the viscodes to be sampled')
    args = parser.parse_args()

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

    main_estimate_dpi(measurements, args.visits, fitter)

    # Get RCDs
    # rcds = data_handler.get_mmse_decline_as_dict()
    # main_estimate_dpi_dpr(measurements, viscodes, fitter, rcds)
    # main_evaluate_scalings(measurements, viscodes, fitter, rcds)


def main_estimate_dpi(measurements, viscodes, fitter):
    # Test all available subjects
    dpis = []
    # progresses = []
    diagnoses = []
    for rid in measurements:
        print log.INFO, 'Estimating DPI for subject {0}...'.format(rid)
        try:
            # progress = measurements[rid]['bl']['progress']
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

        # print log.RESULT, 'Estimated DPI: {0}, Progress: {1}'.format(dpi, progress)
        print log.RESULT, 'Estimated DPI: {0}, Diagnosis: {1}'.format(dpi, diagnosis)
        if dpi is not None:
            dpis.append(dpi)
            # progresses.append(progress)
            diagnoses.append(diagnosis)

    # rms_error = np.sqrt(np.sum(np.square(np.array(progresses) - np.array(dpis))) / len(dpis))
    # mean_error = np.sum(np.abs(np.array(progresses) - np.array(dpis))) / len(dpis)
    # print log.RESULT, 'Mean error: {0}, RMS error: {1}'.format(mean_error, rms_error)

    # Plot the results
    # plot_correlation(dpis, progresses)
    misclassified = 0
    for dpi, diagnosis in zip(dpis, diagnoses):
        if dpi > 0 and diagnosis == 0.5 or dpi <= 0 and diagnosis == 1.0:
            misclassified += 1
    print log.RESULT, 'Accuracy: {0}'.format(1.0 - float(misclassified) / float(len(dpis)))


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


def plot_correlation(dpis, progresses):
    min_progress = np.min(progresses)
    max_progress = np.max(progresses)
    # min_dpi = np.min(progresses)
    # max_dpi = np.max(progresses)
    plt.title('Correlation between progress and estimated DPI')
    plt.xlabel('Estimated DPI')
    plt.ylabel('Disease progression relative to point of conversion')
    # plt.xlim(np.min(dpis), np.max(dpis))
    # plt.ylim(np.min(progresses), np.max(progresses))
    plt.plot([min_progress, max_progress],
             [min_progress, max_progress],
             color='0.5', linestyle='--')
    plt.scatter(dpis, progresses, alpha=0.5)
    plt.show()


def plot_rcds(dpis, dprs, rcds):
    plt.title('DPI and DPR with RCD')
    plt.xlabel('Estimated DPI')
    plt.ylabel('Estimated DPR')
    plt.scatter(dpis, dprs, c=rcds, cmap=aplt.progression_cmap, s=50, linewidths=0)
    plt.show()


if __name__ == '__main__':
    main()
