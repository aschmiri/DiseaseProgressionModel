#! /usr/bin/env python2.7
import numpy as np
import argparse
import matplotlib.pyplot as plt
from common import log as log
from common import adni_tools as adni
from common import adni_plot as aplt
from common import vgam as vgam


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser.add_argument('-i', '--iteration', type=int, default=2, help='the refinement iteration')
    # parser = vgam.add_common_arguments(parser)
    args = parser.parse_args()

    # main_estimate_dpi(args)
    # main_estimate_dpi_dpr(args)
    # main_evaluate_scalings(args)
    main_single_biomarker_ranking(args)


def main_estimate_dpi(args):
    # Define parameters for test
    biomarkers = ['MMSE', 'CDRSB', 'Right Amygdala', 'Left Amygdala', 'Right Lateral Ventricle', 'Left Lateral Ventricle']
    # biomarkers = ['MMSE', 'CDRSB', 'D1', 'D2', 'D3']
    # biomarkers = ['CDRSB']
    # biomarkers = ['Left PCu precuneus']
    # biomarkers = ['Right Amygdala', 'Left Amygdala', 'Right Lateral Ventricle', 'Left Lateral Ventricle', 'Right Hippocampus', 'Left Hippocampus']
    # biomarkers = adni.volume_names
    # biomarkers = adni.cog_score_names
    # viscodes = ['bl']
    viscodes = ['bl', 'm12', 'm24']

    # Collect data for test
    data_files = vgam.get_data_files(args)
    measurements = vgam.get_measurements_as_dict(data_files)

    data_folders = vgam.get_data_folders(args)
    densities = vgam.get_pfds_as_dict(data_folders, biomarkers=biomarkers)

    # Test all available subjects
    dpis = []
    progresses = []
    for rid in measurements:
        print log.INFO, 'Estimating DPI for subject {0}...'.format(rid)
        try:
            progress = measurements[rid]['bl']['progress']
        except Exception:
            continue

        if not set(viscodes).issubset(set(measurements[rid].keys())):
            print log.WARNING, 'Not all viscodes {0} available for subject {1}!'.format(viscodes, rid)
            continue

        samples = {}
        for viscode in viscodes:
            samples.update({viscode: measurements[rid][viscode]})
        dpi = vgam.get_dpi_for_samples(densities, samples, biomarkers=biomarkers)

        print log.RESULT, 'Estimated DPI: {0}, Progress: {1}'.format(dpi, progress)
        dpis.append(dpi)
        progresses.append(progress)

    rms_error = np.sqrt(np.sum(np.square(np.array(progresses) - np.array(dpis))) / len(dpis))
    mean_error = np.sum(np.abs(np.array(progresses) - np.array(dpis))) / len(dpis)
    print log.RESULT, 'Mean error: {0}, RMS error: {1}'.format(mean_error, rms_error)

    # Plot the results
    plot_correlation(dpis, progresses)


def main_estimate_dpi_dpr(args):
    # Define parameters for test
    biomarkers = ['MMSE', 'CDRSB', 'Right Amygdala', 'Left Amygdala', 'Right Lateral Ventricle', 'Left Lateral Ventricle']
    # biomarkers = ['CDRSB']
    # biomarkers = ['Right Amygdala', 'Left Amygdala', 'Right Lateral Ventricle', 'Left Lateral Ventricle', 'Right Hippocampus', 'Left Hippocampus']
    # biomarkers = ['Right Amygdala', 'Left Amygdala', 'Right Lateral Ventricle', 'Left Lateral Ventricle']
    # biomarkers = adni.volume_names
    # biomarkers = adni.cog_score_names
    # viscodes = ['bl']
    viscodes = ['bl', 'm12', 'm24']

    # Collect data for test
    data_files = vgam.get_data_files(args)
    measurements = vgam.get_measurements_as_dict(data_files)

    data_folders = vgam.get_data_folders(args)
    densities = vgam.get_pfds_as_dict(data_folders, biomarkers=biomarkers)

    rcds = vgam.get_rcd_as_dict(measurements)

    # Test all available subjects
    dpis = []
    dprs = []
    progresses = []
    rcdnum = []
    for rid in measurements:
        print log.INFO, 'Estimating DPI and DPR for subject {0}...'.format(rid)
        try:
            progress = measurements[rid]['bl']['progress']
        except Exception:
            continue

        if not set(viscodes).issubset(set(measurements[rid].keys())):
            print log.WARNING, 'Not all viscodes {0} available for subject {1}!'.format(viscodes, rid)
            continue

        samples = {}
        for viscode in viscodes:
            samples.update({viscode: measurements[rid][viscode]})
        dpi, dpr = vgam.get_dpi_dpr_for_samples(densities, samples, biomarkers=biomarkers)

        print log.RESULT, 'Estimated DPI: {0}, DPR: {1}, Progress: {2}'.format(dpi, dpr, progress)
        dpis.append(dpi)
        dprs.append(dpr)
        progresses.append(progress)
        # rcdnum.append(1 if rcds[rid] else 0)
        rcdnum.append(rcds[rid])

    rms_error = np.sqrt(np.sum(np.square(np.array(progresses) - np.array(dpis))) / len(dpis))
    mean_error = np.sum(np.abs(np.array(progresses) - np.array(dpis))) / len(dpis)
    print log.RESULT, 'Mean error: {0}, RMS error: {1}'.format(mean_error, rms_error)

    # Plot the results
    plot_rcds(dpis, dprs, rcdnum)


def main_evaluate_scalings(args):
    biomarkers = ['MMSE']

    # Collect data for test
    data_files = vgam.get_data_files(args)
    measurements = vgam.get_measurements_as_dict(data_files)

    input_folders = vgam.get_data_folders(args)
    measurements = vgam.get_scaled_measurements(input_folders, measurements, biomarkers=biomarkers)

    rcds = vgam.get_rcd_as_dict(measurements)

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


def main_single_biomarker_ranking(args):
    # Collect data for test
    data_files = vgam.get_data_files(args)
    measurements = vgam.get_measurements_as_dict(data_files)
    data_folders = vgam.get_data_folders(args)

    # Compute error for each biomarker
    mean_errors = []
    for biomarker in adni.manifold_coordinate_names:
        densities = vgam.get_pfds_as_dict(data_folders, biomarkers=[biomarker])

        # Test all available subjects
        print log.INFO, 'Computing RMS for {0}...'.format(biomarker)
        rms_error = 0
        num_errors = 0
        for rid in measurements:
            for viscode in measurements[rid]:
                progress = measurements[rid][viscode]['progress']

                samples = {}
                samples.update({viscode: measurements[rid][viscode]})
                dpi = vgam.get_dpi_for_samples(densities, samples, biomarkers=[biomarker])
                rms_error += np.square(dpi - progress)
                num_errors += 1

        rms_error /= num_errors
        rms_error = np.sqrt(rms_error)
        mean_errors.append(rms_error)
        print log.RESULT, 'RMS error for {0}: {1} ({2})'.format(biomarker, rms_error, num_errors)

    # Sort results according to error
    mean_errors = np.array(mean_errors)
    biomarker_names = np.array(adni.biomarker_names)

    args = np.argsort(mean_errors)
    mean_errors = mean_errors[args]
    biomarker_names = biomarker_names[args]

    for name, error in zip(biomarker_names, mean_errors):
        print name, ':', error


def plot_correlation(dpis, progresses):
    min_progress = np.min(progresses)
    max_progress = np.max(progresses)

    plt.title('Correlation between progress and estimated DPI')
    plt.xlabel('Estimated DPI')
    plt.ylabel('Disease progression relative to point of conversion')
    plt.xlim(min_progress, 20)
    plt.ylim(min_progress, 5)
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
