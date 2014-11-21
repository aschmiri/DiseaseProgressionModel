#! /usr/bin/env python2.7
import os.path
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from common import log as log
from common import plotting_tools as pt
from common.datahandler import DataHandler
from common.progressionmodel import ProgressionModel


def main():
    parser = argparse.ArgumentParser(description='Estimate model curves for biomarkers using VGAM.')
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-e', '--extrapolator', type=str, choices=['lin', 'sqrt', 'exp'], default='exp', help='the type of extrapolator')
    parser.add_argument('--recompute_errors', action='store_true', help='recompute the matrix containing the fitting errors')
    parser.add_argument('--search_range', nargs=3, default=(1000, 5000, 10), help='the range in which the offset is sought')
    args = parser.parse_args()

    # Get the data files and biomarkers
    data_handler_joint = DataHandler.get_data_handler(method=args.method,
                                                      biomarkers=args.biomarkers,
                                                      phase='joint')
    biomarkers, offsets, errors, descriminativeness, overlap = get_fitting_data(args, data_handler_joint)

    # Plot single biomarker fits
    fig, ax = plt.subplots()
    pt.setup_axes(plt, ax, xgrid=False)
    ax.set_title('Optimal offset between CN/MCI and MCI/AD models')
    ax.set_xlabel('Offset (days)')
    ax.set_ylabel('Fitting error')
    for i, biomarker in enumerate(biomarkers):
        if descriminativeness[i] > 0.4:
            print log.RESULT, 'Min error for {0} at {1}'.format(biomarker, offsets[np.argmin(errors[i, :])])
            ax.plot(offsets, errors[i, :], label=biomarker, linestyle='--')

    # Get optimal offset
    mean_errors = np.mean(errors, 0)
    weighted_mean_errors = np.dot(errors.T, descriminativeness) / np.sum(descriminativeness)

    # Plot joint fit
    ax.plot(offsets, mean_errors, label='Mean', linewidth=2, color='g')
    ax.plot(offsets, weighted_mean_errors, label='Weighted mean', linewidth=2, color='r')

    # Get and lot optimal offset
    optimal_offset = offsets[np.argmin(mean_errors)]
    optimal_offset_weighted = offsets[np.argmin(weighted_mean_errors)]
    print log.RESULT, 'Optimal threshold: {0}'.format(optimal_offset)
    print log.RESULT, 'Optimal threshold (weighted): {0}'.format(optimal_offset_weighted)
    ax.axvline(optimal_offset, linestyle=':', color='g')
    ax.axvline(optimal_offset_weighted, linestyle=':', color='r')

    # Plot overlap
    ax.axvline(overlap, color='0.15', linestyle=':')

    ax.legend()
    plt.show()
    plt.close(fig)


def get_fitting_data(args, data_handler_joint):
    biomarkers = data_handler_joint.get_biomarker_names()
    offsets = range(args.search_range[0], args.search_range[1], args.search_range[2])
    errors_file = os.path.join(data_handler_joint.get_eval_folder(),
                               'offset_errors_{0}.p'.format(args.extrapolator))
    if os.path.isfile(errors_file) and not args.recompute_errors:
        print log.INFO, 'Reading errors estimations from file {0}...'.format(errors_file)
        (errors, descriminativeness, overlap) = pickle.load(open(errors_file, 'rb'))
    else:
        data_handler_1 = DataHandler.get_data_handler(method=args.method,
                                                      biomarkers=args.biomarkers,
                                                      phase='cnmci')
        data_handler_2 = DataHandler.get_data_handler(method=args.method,
                                                      biomarkers=args.biomarkers,
                                                      phase='mciad')

        errors = np.zeros((len(biomarkers), len(offsets)))
        descriminativeness = np.zeros(len(biomarkers))
        overlap = []
        for i, biomarker in enumerate(biomarkers):
            # Get error matrix for all biomarkers and offsets
            model_file_1 = data_handler_1.get_model_file(biomarker)
            model_file_2 = data_handler_2.get_model_file(biomarker)
            if os.path.isfile(model_file_1) and os.path.isfile(model_file_2):
                print log.INFO, 'Analysing {0}...'.format(biomarker)

                # Get discriminativeness for all biomarkers as a scaling factor
                eval_file_1 = model_file_1.replace('.csv', '_eval_cover.csv')
                eval_file_2 = model_file_2.replace('.csv', '_eval_cover.csv')
                if os.path.isfile(eval_file_1) and os.path.isfile(eval_file_2):
                    descriminate_1 = np.mean(mlab.csv2rec(eval_file_1)['error'])
                    descriminate_2 = np.mean(mlab.csv2rec(eval_file_2)['error'])
                    descriminativeness[i] = 0.5 * (descriminate_1 + descriminate_2)
                else:
                    print log.WARNING, 'Evaluation file missing for {0}'.format(biomarker)
                    continue

                # Initialise models
                model_1 = ProgressionModel(biomarker, model_file_1, extrapolator=args.extrapolator)
                model_2 = ProgressionModel(biomarker, model_file_2, extrapolator=args.extrapolator)

                # Assemble errors for each offset
                min_val_1, max_val_1 = model_1.get_value_range([0.1, 0.9])
                min_val_2, max_val_2 = model_2.get_value_range([0.1, 0.9])
                values = np.linspace(min(min_val_1, min_val_2), max(max_val_1, max_val_2), 250)
                values_delta = (values.max() - values.min()) / len(values)
                for j, offset in enumerate(offsets):
                    dens_11 = np.array(model_1.get_density_distribution(values, offset + model_2.min_progress))
                    dens_12 = np.array(model_2.get_density_distribution(values, model_2.min_progress))

                    dens_21 = np.array(model_1.get_density_distribution(values, model_1.max_progress))
                    dens_22 = np.array(model_2.get_density_distribution(values, -offset + model_1.max_progress))

                    errors[i, j] = 0.5 * values_delta * (np.sum(np.abs(dens_11 - dens_12)) + np.sum(np.abs(dens_21 - dens_22)))

                # Get overlap
                overlap.append(model_1.max_progress - model_2.min_progress)

        overlap = np.mean(overlap)
        print log.INFO, 'Saving errors to file {0}...'.format(errors_file)
        pickle.dump((errors, descriminativeness, overlap), open(errors_file, 'wb'))

    return biomarkers, offsets, errors, descriminativeness, overlap

if __name__ == '__main__':
    main()
