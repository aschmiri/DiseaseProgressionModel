#! /usr/bin/env python2.7
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from common import log as log
from common import plotting_tools as pt
from common.progressionmodel import ProgressionModel
from common.synthmodel import SynthModel


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--biomarker', default='cdrsb', choices=['cdrsb', 'mmse', 'hipp'],
                        help='name of the biomarker to be evaluated')
    parser.add_argument('--num_samples', type=int, default=200, help='number of samples')
    parser.add_argument('--show_plots', action='store_true', default=False,
                        help='show the plots with a default filename')
    args = parser.parse_args()

    # evaluate_curves('synth_cdrsb', num_samples=args.num_samples, show_plots=args.show_plots)

    # for num_samples in range(100, 2000, 100):
    #     evaluate_curves('synth_cdrsb', num_samples=num_samples, show_plots=args.show_plots)

    for csig in [0, 200, 400]:
        evaluate_curves(args.biomarker, num_samples=1000, show_plots=args.show_plots, csig=csig)


def evaluate_curves(biomarker_name, num_samples=200, show_plots=False, csig=0):
    biomarker = 'synth_{0}'.format(biomarker_name)
    print log.INFO, 'Evaluating {0} for {1} samples, csig={2}...'.format(biomarker, num_samples, csig)
    donohue_model_path = '/Development/DiseaseProgressionModel/models/donohue/'
    vgam_model_path = '/Development/DiseaseProgressionModel/models/synth/'

    # Setup plot
    if show_plots:
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
        pt.setup_axes(plt, ax, xgrid=False, ygrid=False)
        ax.set_title('')
        ax.set_xlabel('')
    else:
        fig = None
        ax = None

    # Initialise values
    offset_donohue = 0  # 182.5
    errors_donohue = []
    errors_vgam_mean = []
    errors_vgam_median = []

    # Get real curve values
    progress_linspace = np.linspace(-1500, 1500)
    mean_curve = [SynthModel.get_mean_value(biomarker, p) for p in progress_linspace]
    median_curve = [SynthModel.get_distributed_value(biomarker, p, cdf=0.5) for p in progress_linspace]

    # Plot synthetic model curve
    if show_plots:
        progress_linspace_synth = np.linspace(-2500, 2500, 100)
        quantiles = [0.1, 0.25, 0.5, 0.75, 0.9]
        alphas = [0.4, 0.7, 1.0, 0.7, 0.4]
        for quantile, alpha in zip(quantiles, alphas):
            curve_synth = [SynthModel.get_distributed_value(biomarker, p, cdf=quantile)
                           for p in progress_linspace_synth]
            ax.plot(progress_linspace_synth, curve_synth, color='b', alpha=alpha)
        curve_synth = [SynthModel.get_mean_value(biomarker, p) for p in progress_linspace_synth]
        ax.plot(progress_linspace_synth, curve_synth, '--', color='b')

    # Get values for mean calculation
    end_values = [SynthModel.get_distributed_value(biomarker, progress_linspace[0], cdf=0.001),
                  SynthModel.get_distributed_value(biomarker, progress_linspace[-1], cdf=0.001),
                  SynthModel.get_distributed_value(biomarker, progress_linspace[0], cdf=0.999),
                  SynthModel.get_distributed_value(biomarker, progress_linspace[-1], cdf=0.999)]
    values = np.linspace(min(end_values), max(end_values), 100)
    for run in range(100):
        # Get Donohue model
        donohue_file = os.path.join(donohue_model_path,
                                    'population_value-{0}_csig{1}_run{2}.csv'.format(biomarker_name, csig, run))

        r = mlab.csv2rec(donohue_file)
        progrs = r[r.dtype.names[0]] - offset_donohue
        vals = r[r.dtype.names[1]]

        curve_donohue = []
        progr_donohue = []
        for p in progress_linspace:
            if progrs[0] < p < progrs[-1]:
                i = 1
                while p > progrs[i]:
                    i += 1
                progr_donohue.append(float(progrs[i]))
                curve_donohue.append(float(vals[i]))
            else:
                print log.WARNING, 'Model scope too small... skipping!'
                continue

        # Get VGAM model
        if csig == 0:
            vgam_model_file = os.path.join(vgam_model_path,
                                           '{0}_model_{1}_longitudinal_{2}.csv'.format(biomarker,
                                                                                       num_samples, run))
        else:
            vgam_model_file = os.path.join(vgam_model_path,
                                           '{0}_model_{1}_longitudinal_csig{2}.0_{3}.csv'.format(biomarker,
                                                                                                 num_samples,
                                                                                                 csig, run))

        pm = ProgressionModel(biomarker, vgam_model_file)
        curve_vgam_median = pm.get_quantile_curve(progress_linspace, 0.5)
        curve_vgam_mean = [np.sum(pm.get_density_distribution(values, p) * values /
                                  np.sum(pm.get_density_distribution(values, p))) for p in progress_linspace]

        # Calculate errors
        errors_donohue.append(np.mean(np.abs(np.array(curve_donohue) - np.array(mean_curve))))
        errors_vgam_mean.append(np.mean(np.abs(np.array(curve_vgam_mean) - np.array(mean_curve))))
        errors_vgam_median.append(np.mean(np.abs(np.array(curve_vgam_median) - np.array(median_curve))))

        if show_plots:
            ax.plot(progr_donohue, curve_donohue, '--', color='g', alpha=0.2, linewidth=2)
            # ax.plot(progress_linspace, curve_vgam_median, '-', color='r', alpha=0.2, linewidth=2)
            ax.plot(progress_linspace, curve_vgam_mean, '--', color='r', alpha=0.2, linewidth=2)

    print log.RESULT, 'Donohue (mean):', np.mean(errors_donohue), np.var(errors_donohue)
    print log.RESULT, 'VGAM    (mean):', np.mean(errors_vgam_mean), np.var(errors_vgam_mean)
    print log.RESULT, 'VGAM  (median):', np.mean(errors_vgam_median), np.var(errors_vgam_median)

    # Draw or save the plot
    if show_plots:
        plt.tight_layout()
        plt.show()
        plt.close(fig)


if __name__ == '__main__':
    main()
