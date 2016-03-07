#! /usr/bin/env python2.7
# coding=utf-8
"""
A meant for the VPH-Dare@IT platform script to
- read a parameter file containing all the biomarker values of a subject
- fit these biomarkers to the (already generated) models to estimate the subject’s disease progression and
- plot graphs of the respective models with a line indicating the estimated progression


:author:     Alexander Schmidt-Richberg
:copyright:  2016 Imperial College London. All rights reserved.
:contact:    a.schmidt-richberg@imperial.ac.uk
"""
import argparse
import os
import csv
import numpy as np

import matplotlib.pyplot as plt

from common import log as log
from common import plotting_tools as pt
from common.datahandler import DataHandler
from common.modelfitter import ModelFitter
from common.progressionmodel import MultiBiomarkerProgressionModel, ProgressionModel
from common.evaluation_tools import estimate_dpis, estimate_dpis_dprs


def read_measurements_from_cvs(filename):
    """
    Created a dict from the sample measurements file. For compatibility with the library, the dict has to have the
    { <rid> : { <viscode> : { DX.scan : <diagnosis> }
                            { scantime : <days after bl> }
                            { <biomarker1> : <volume> }
                        ... }
              { <viscode> : ... }}

    :param filename: filename of the *.csv file
    :rtype: dict
    :return: the generated dict with the measurements
    """
    scantime_dict = {'bl': 0, 'm12': 365, 'm24': 730, 'm36': 1095}

    biomarkers = set()
    measurements = {0: {}}
    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        visits = reader.fieldnames[1:]
        for visit in visits:
            measurements[0].update({visit: {'scantime': scantime_dict[visit], 'DX.scan': 'UNKNOWN'}})

        for row in reader:
            biomarker = row['Biomarker Name']
            if biomarker in DataHandler.get_all_biomarker_names():
                for visit in visits:
                    measurements[0][visit].update({biomarker: float(row[visit])})
                    biomarkers.add(biomarker)

    return measurements, list(biomarkers)


def plot_biomarker(data_handler, biomarker, measurements, dpi, dpr):
    """
    Plot the model of one biomarker with the fitted values

    :param data_handler: the data handler
    :param biomarker: the biomarker to plot
    :param measurements: the measurements containing the biomarker samples of one subject
    :param dpi: the estimated DPI
    :param dpr: the estimated DPR
    """
    model_file = data_handler.get_model_file(biomarker)
    if not os.path.isfile(model_file):
        print log.ERROR, 'Model file not found: {0}'.format(model_file)
        return

    print log.INFO, 'Generating plot for {0}...'.format(biomarker)

    #
    # Read model
    #
    pm = ProgressionModel(biomarker, model_file)
    progress_extrapolate = 0.3 * (pm.max_progress - pm.min_progress)
    min_progress_extrapolate = int(pm.min_progress - progress_extrapolate)
    max_progress_extrapolate = int(pm.max_progress + progress_extrapolate)
    progress_linspace_ex1 = np.linspace(min_progress_extrapolate, pm.min_progress, 20)
    progress_linspace_int = np.linspace(pm.min_progress, pm.max_progress, 60)
    progress_linspace_ex2 = np.linspace(pm.max_progress, max_progress_extrapolate, 20)

    #
    # Setup plot
    #
    biomarker_string = pt.get_biomarker_string(biomarker)
    figure_width = 6
    fig = plt.figure(figsize=(figure_width, 5))
    ax1 = plt.subplot(1, 1, 1)
    pt.setup_axes(plt, ax1, xgrid=False, ygrid=False)
    ax1.set_title('Model for {0} with fitted sample values'.format(biomarker_string))
    ax1.set_xlabel('Disease progress (days before/after conversion to MCI)')
    ax1.set_ylabel(DataHandler.get_biomarker_unit(biomarker))
    ax1.set_xlim(min_progress_extrapolate, max_progress_extrapolate)

    #
    # Plot the percentile curves of the fitted model
    #
    ax1.axvline(pm.min_progress, color='0.15', linestyle=':')
    ax1.axvline(pm.max_progress, color='0.15', linestyle=':')

    quantiles = [0.1, 0.25, 0.5, 0.75, 0.9]
    grey_values = ['0.4', '0.2', '0', '0.2', '0.4']
    for grey_value, quantile in zip(grey_values, quantiles):
        curve_int = pm.get_quantile_curve(progress_linspace_int, quantile)
        ax1.plot(progress_linspace_int, curve_int, color=grey_value)

        curve_ex1 = pm.get_quantile_curve(progress_linspace_ex1, quantile)
        curve_ex2 = pm.get_quantile_curve(progress_linspace_ex2, quantile)
        ax1.plot(progress_linspace_ex1, curve_ex1, '--', color=grey_value)
        ax1.plot(progress_linspace_ex2, curve_ex2, '--', color=grey_value)

        label = 'q = {0}'.format(quantile * 100)
        ax1.text(progress_linspace_int[-1] + 100, curve_int[-1], label, fontsize=10)

    #
    # Plot points
    #
    progr_points = []
    value_points = []
    diagn_points = []
    for visit in measurements[0]:
        progress = measurements[0][visit]['scantime'] * dpr + dpi
        value = measurements[0][visit][biomarker]
        progr_points.append(progress)
        value_points.append(value)
        diagn_points.append(1.0)
        ax1.axvline(progress, color='b', linestyle='--')
        ax1.text(progress + 150, value, visit, color='b', fontsize=10)

    ax1.scatter(progr_points, value_points, s=25.0, color='b', edgecolor='none',
                vmin=0.0, vmax=1.0, alpha=0.9)

    #
    # Draw or save the plot
    #
    plt.tight_layout()
    plt.show()
    plt.close(fig)


def main():
    # Parse input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--estimate_dpr', action='store_true', help='recompute the dpis estimations')
    parser.add_argument('--samples_file', type=str, default='measurements_sample.csv',
                        help='recompute the dpis estimations')
    args = parser.parse_args()

    # Read the measurements as dict from the csv file
    measurements, biomarkers = read_measurements_from_cvs(args.samples_file)
    visits = measurements[0].keys()

    # Get estimates
    data_handler = DataHandler.get_data_handler(method='all',
                                                biomarkers=biomarkers,
                                                phase='joint')

    # Setup model
    model = MultiBiomarkerProgressionModel()
    for biomarker in biomarkers:
        model_file = data_handler.get_model_file(biomarker)
        model.add_model(biomarker, model_file)
    fitter = ModelFitter(model)

    # Estimate dpis (and dprs) and save data
    if args.estimate_dpr:
        rids, diagnoses, dpis, dprs = estimate_dpis_dprs(measurements, visits, fitter, phase='joint')
    else:
        rids, diagnoses, dpis = estimate_dpis(measurements, visits, fitter, phase='joint')
        dprs = np.ones(len(dpis)).tolist()

    # Plot the models with the fitted samples
    for biomarker in biomarkers:
        plot_biomarker(data_handler, biomarker, measurements, dpis[0], dprs[0])


if __name__ == '__main__':
    main()
