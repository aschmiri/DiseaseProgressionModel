#! /usr/bin/env python
# print __doc__
import os.path
import math
import numpy as np
import matplotlib.pyplot as plt
from src.common import adni_tools as adni
from src.common import vgam as vgam


def main():
    # Define parameters for test
    # biomarkers = ['MMSE', 'CDRSB', 'Right Amygdala', 'Left Amygdala', 'Right Lateral Ventricle', 'Left Lateral Ventricle']
    biomarkers = ['Left PCu precuneus']
    # biomarkers = ['Right Amygdala', 'Left Amygdala', 'Right Lateral Ventricle', 'Left Lateral Ventricle', 'Right Hippocampus', 'Left Hippocampus']
    # biomarkers = adni.volume_names
    # biomarkers = adni.cog_score_names
    viscodes = ['bl']
    # viscodes = ['bl', 'm12', 'm24']

    # Collect data for test
    data_file = os.path.join(adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv')
    measurements = vgam.get_measurements_as_collection(data_file)
    densities = vgam.get_pfds_as_collection(biomarkers=biomarkers)

    # Test all available subjects
    dpis = []
    progresses = []
    mean_error = 0
    for rid in measurements:
        print 'Estimating DPI for subject', rid
        progress = measurements[rid]['bl']['progress']

        samples = {}
        for viscode in viscodes:
            samples.update({viscode: measurements[rid][viscode]})
        dpi = evaluate_samples(densities, samples, biomarkers=biomarkers)

        print 'Estimated DPI', dpi, '; Progress:', progress
        dpis.append(dpi)
        progresses.append(progress)
        mean_error += np.abs(dpi - progress)

    mean_error /= len(measurements)
    print 'Mean error:', mean_error

    plot_correlation(dpis, progresses)


def main_single_biomarker_ranking():
    # Collect data for test
    data_file = os.path.join(adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv')
    measurements = vgam.get_measurements_as_collection(data_file)

    # Compute error for each biomarker
    mean_errors = []
    for biomarker in adni.biomarker_names:
        densities = vgam.get_pfds_as_collection(biomarkers=[biomarker])

        # Test all available subjects
        mean_error = 0
        num_errors = 0
        for rid in measurements:
            for viscode in measurements[rid]:
                progress = measurements[rid][viscode]['progress']

                samples = {}
                samples.update({viscode: measurements[rid][viscode]})
                dpi = evaluate_samples(densities, samples, biomarkers=[biomarker])
                mean_error += np.square(dpi - progress)
                num_errors += 1

        mean_error /= num_errors
        mean_error = math.sqrt(mean_error)
        mean_errors.append(mean_error)
        print 'Mean error for {0}: {1} ({2})'.format(biomarker, mean_error, num_errors)

    # Sort results according to error
    mean_errors = np.array(mean_errors)
    biomarker_names = np.array(adni.biomarker_names)

    args = np.argsort(mean_errors)
    mean_errors = mean_errors[args]
    biomarker_names = biomarker_names[args]

    for name, error in zip(biomarker_names, mean_errors):
        print name, ':', error


def evaluate_sample_at_progress(densities, sample, progress, biomarkers=adni.biomarker_names):
    prob_sample = 1.0

    prog_prev = math.floor(progress)
    prog_next = math.ceil(progress)
    prog_offset = progress - prog_prev

    for biomarker in biomarkers:
        if biomarker in densities:
            values = densities[biomarker]['values']
            value_sample = sample[biomarker]

            if value_sample == None:
                print 'WARNING: sample has no value for', biomarker
            else:
                # Find value in probability list
                i = 0
                while values[i] < value_sample and i < len(values):
                    i += 1
                i_prev = i - 1 if i > 0 else 0
                i_next = i if i < len(values) else len(values)

                # Get factor for value
                value_prev = values[i_prev]
                value_next = values[i_next]
                factor = (value_sample - value_prev) / (value_next - value_prev)

                # Interpolate probability
                prob_prev = densities[biomarker][prog_prev][i_prev]
                prob_next = densities[biomarker][prog_prev][i_next]
                prob_sample_prev = factor * prob_next + (1 - factor) * prob_prev

                if prog_offset == 0.0:
                    prob_sample *= prob_sample_prev
                else:
                    prob_prev = densities[biomarker][prog_next][i_prev]
                    prob_next = densities[biomarker][prog_next][i_next]
                    prob_sample_next = factor * prob_next + (1 - factor) * prob_prev

                    # Interpolate between
                    prob_sample *= (1 - prog_offset) * prob_sample_prev + prog_offset * prob_sample_next

    return prob_sample


def evaluate_samples(densities, samples, biomarkers=adni.biomarker_names):
    max_scantime = np.max([samples[viscode]['scantime'] for viscode in samples])
    progress_steps = np.arange(vgam.MIN_PROGRESS, vgam.MAX_PROGRESS - max_scantime, 0.5)
    probs = []
    for progress in progress_steps:
        prob = 1.0
        for viscode in samples:
            offset = samples[viscode]['scantime']
            prob *= evaluate_sample_at_progress(densities, samples[viscode], progress + offset, biomarkers=biomarkers)
        probs.append(prob)
    dpi = progress_steps[np.argmax(probs)]
    return dpi


def plot_correlation(dpis, progresses):
    plt.title('Correlation between progress and estimated DPI')
    plt.xlabel('Disease progression relative to point of conversion')
    plt.ylabel('Estimated DPI')
    plt.xlim(vgam.MIN_PROGRESS, 20)
    plt.ylim(vgam.MIN_PROGRESS, 5)
    plt.plot([vgam.MIN_PROGRESS, vgam.MAX_PROGRESS],
             [vgam.MIN_PROGRESS, vgam.MAX_PROGRESS],
             color='0.5', linestyle='--')
    plt.scatter(dpis, progresses, alpha=0.5)
    plt.show()


if __name__ == '__main__':
    main_single_biomarker_ranking()
    # main()
