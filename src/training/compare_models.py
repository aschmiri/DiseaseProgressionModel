#! /usr/bin/env python2.7
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from common import log as log
from common.progressionmodel import ProgressionModel
from common.datahandler import DataHandler


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('-p', '--phase', default='mciad', choices=DataHandler.get_phase_choices(), help='the phase for which the model is to be trained')
    parser.add_argument('-e', '--extrapolator', type=str, choices=['lin', 'sqrt', 'exp'], default='exp', help='the type of extrapolator')
    args = parser.parse_args()

    data_handler = DataHandler.get_data_handler(method=args.method,
                                                biomarkers=args.biomarkers,
                                                phase=args.phase)

    biomarkers = data_handler.get_biomarker_names()
    if args.method == 'joint':
        offsets = np.linspace(500, 3000, 26)
    else:
        offsets = np.linspace(-1000, 1000, 21)
    all_diffs = np.zeros((len(offsets), len(biomarkers)))

    for i, biomarker in enumerate(biomarkers):
        diffs = get_model_differences(args, data_handler, biomarker, offsets)
        all_diffs[:, i] = diffs
        print biomarker, offsets[np.argmin(diffs)]

    optimum_index = np.argmin(np.mean(all_diffs, axis=1))

    print 'all', offsets[optimum_index]

    mins = all_diffs[optimum_index, :]  # np.min(all_diffs, axis=0)
    indices = np.argsort(mins)
    for i in indices:
        print biomarkers[i], mins[i]

    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)
    ax1.plot(offsets, all_diffs, color='r')
    ax1.plot(offsets, np.mean(all_diffs, axis=1), color='b')
    plt.show()
    plt.close(fig)

def get_model_differences(args, data_handler, biomarker, offsets):
    print log.INFO, 'Comparing models for {0}...'.format(biomarker)

    model_file = data_handler.get_model_file(biomarker)
    print model_file
    if not os.path.isfile(model_file):
        print log.ERROR, 'Model file not found: {0}'.format(model_file)
        return
    donohue_model_file = os.path.join(data_handler._conf.models_folder,
                                      'denohue', 'population_{0}.csv'.format(biomarker.replace(' ', '.')))
    if not os.path.isfile(donohue_model_file):
        print log.ERROR, 'Donohue Model file not found: {0}'.format(model_file)
        return

    # Read Donohue model
    r = mlab.csv2rec(donohue_model_file)
    progrs = r[r.dtype.names[0]] * 30.44
    progrs = progrs[::100]
    vals_donohue = r[r.dtype.names[1]]
    vals_donohue = vals_donohue[::100]

    # Read my model
    pm = ProgressionModel(biomarker, model_file, extrapolator=args.extrapolator)
    diffs = np.empty(len(offsets))
    for i, offset in enumerate(offsets):
        vals_mine = pm.get_quantile_curve(progrs + offset, 0.5)
        normalizer = max(np.max(vals_mine), np.max(vals_donohue))
        diffs[i] = np.mean(np.abs(vals_donohue - vals_mine)) /normalizer

    return diffs

if __name__ == '__main__':
    main()
