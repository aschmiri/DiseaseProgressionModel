#! /usr/bin/env python
# print __doc__

import os.path
import math
import numpy as np
import matplotlib.pyplot as plt
from src.common import adni_tools as adni
from src.common import vgam as vgam


def main():
    biomarkers = ['MMSE', 'CDRSB', 'Right Amygdala', 'Left Amygdala', 'Right Lateral Ventricle', 'Left Lateral Ventricle']
    #biomarkers = ['MMSE', 'CDRSB']
    #biomarkers = ['Right Amygdala', 'Left Amygdala', 'Right Lateral Ventricle', 'Left Lateral Ventricle', 'Right Hippocampus', 'Left Hippocampus']
    #biomarkers = adni.volume_names
    #biomarkers = adni.cog_score_names

    data_file = os.path.join(adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv')
    measurements = vgam.get_measurements_as_collection(data_file)
    densities = vgam.get_densities_as_collection(biomarkers=biomarkers)

    dpis = []
    progresses = []
    rms = 0
    for rid in measurements:
        print 'Estimating DPI for subject', rid
        progress = measurements[rid]['bl']['progress']

        #sample = measurements[rid]['bl']
        #dpi = evaluate_sample(densities, sample, biomarkers=biomarkers)

        samples = {}
        samples.update({'bl': measurements[rid]['bl']})
        samples.update({'m12': measurements[rid]['m12']})
        dpi = evaluate_samples(densities, samples, biomarkers=biomarkers)

        print 'Estimated DPI', dpi, '; Progress:', progress
        dpis.append(dpi)
        progresses.append(progress)
        rms += np.square(dpi - progress)
    rms = math.sqrt(rms)
    print 'RMS:', rms
    plot_correlation(dpis, progresses)


def plot_correlation(dpis, progresses):
    plt.xlim(vgam.MIN_PROGRESS, 20)
    plt.plot([vgam.MIN_PROGRESS, vgam.MAX_PROGRESS], [vgam.MIN_PROGRESS, vgam.MAX_PROGRESS], color='0.5', linestyle='--')
    plt.scatter(dpis, progresses, alpha=0.5)
    plt.show()


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
                    print biomarker, prog_next, i_prev
                    prob_prev = densities[biomarker][prog_next][i_prev]
                    prob_next = densities[biomarker][prog_next][i_next]
                    prob_sample_next = factor * prob_next + (1 - factor) * prob_prev

                    # Interpolate between
                    prob_sample *= (1 - prog_offset) * prob_sample_prev + prog_offset * prob_sample_next

    return prob_sample


def evaluate_sample(densities, sample, biomarkers=adni.biomarker_names):
    progress_steps = np.arange(vgam.MIN_PROGRESS, vgam.MAX_PROGRESS, 0.5)
    probs = []
    for progress in progress_steps:
        prob = evaluate_sample_at_progress(densities, sample, progress, biomarkers=biomarkers)
        probs.append(prob)

    #plt.plot(progress_steps, probs)
    #plt.show()

    dpi = progress_steps[np.argmax(probs)]
    return dpi


def evaluate_samples(densities, samples, biomarkers=adni.biomarker_names):
    progress_steps = np.arange(vgam.MIN_PROGRESS, vgam.MAX_PROGRESS - 12, 0.5)
    probs = []
    for progress in progress_steps:
        prob = 1.0
        for viscode in samples:
            offset = samples[viscode]['scantime']
            prob *= evaluate_sample_at_progress(densities, samples[viscode], progress + offset, biomarkers=biomarkers)
        probs.append(prob)
    dpi = progress_steps[np.argmax(probs)]
    return dpi

if __name__ == '__main__':
    main()
