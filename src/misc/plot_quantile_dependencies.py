#! /usr/bin/env python2.7
import os.path
import argparse
import numpy as np
import pickle
from scipy.stats import norm
import matplotlib.pyplot as plt
from common import log as log
from common import plotting_tools as pt
from common.datahandler import DataHandler
from common.progressionmodel import ProgressionModel


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--method', choices=DataHandler.get_method_choices(), default='all', help='the method to collect data for')
    parser.add_argument('-b', '--biomarkers', nargs='+', default=None, help='name of the biomarker to be plotted')
    parser.add_argument('--save_plots', action='store_true', default=False, help='save the plots with a default filename')
    args = parser.parse_args()

    # Collect data for test
    data_handler = DataHandler.get_data_handler(method=args.method, biomarkers=args.biomarkers)
    biomarkers = data_handler.get_biomarker_names()
    measurements = data_handler.get_measurements_as_dict(visits=['bl', 'm12'],
                                                         biomarkers=biomarkers,
                                                         select_training_set=True,
                                                         select_complete=True)

    # Setup plotting folder
    eval_folder = DataHandler.make_dir(DataHandler.get_eval_folder(), 'quants')

    # Process all biomarkers
    for biomarker in biomarkers:
        print log.INFO, 'Generating quantile correlation plot for {0}...'.format(biomarker)
        model_file = DataHandler.get_model_file(biomarker)
        pm = ProgressionModel(biomarker, model_file)

        q_file = os.path.join(eval_folder, '{0}.p'.format(biomarker))

        if os.path.isfile(q_file):
            (q_bl, q_m12) = pickle.load(open(q_file, 'rb'))
        else:
            q_bl = []
            q_m12 = []

            for rid in measurements:
                val_bl = measurements[rid]['bl'][biomarker]
                val_m12 = measurements[rid]['m12'][biomarker]

                p_bl = measurements[rid]['bl']['progress']
                p_m12 = measurements[rid]['m12']['progress']

                q_bl.append(pm.approximate_quantile(p_bl, val_bl))
                q_m12.append(pm.approximate_quantile(p_m12, val_m12))

            pickle.dump((q_bl, q_m12), open(q_file, 'wb'))

        # Setup plot
        fig, axs = plt.subplots(1, 2)
        plt.suptitle('Correlation between bl and m12 quantiles')

        # Plot 1
        ax = axs[0]
        pt.setup_axes(plt, ax, yspine=True)
        ax.set_xlabel('Quantile bl')
        ax.set_ylabel('Quantile m12')

        ax.scatter(q_bl, q_m12, edgecolor='none', s=25.0, alpha=0.5)

        # Plot 2
        q_bl = np.array(q_bl)
        q_m12 = np.array(q_m12)

        errors = q_bl - q_m12
        loc, scale = norm.fit(errors, floc=0.0)

        ax = axs[1]
        pt.setup_axes(plt, ax)
        ax.set_xlabel('Difference bl to m12')
        ax.set_ylabel('Probability')
        ax.set_xlim(-1.05, 1.05)
        ax.hist(errors, bins=15, normed=True, histtype='stepfilled', alpha=0.3)
        x = np.linspace(-1.0, 1.0, 100)
        ax.plot(x, norm.pdf(x, loc=loc, scale=scale), color='k')

        # Draw or save the plot
        plt.tight_layout()
        if args.save_plots:
            plot_file = os.path.join(eval_folder, '{0}.pdf'.format(biomarker))
            plt.savefig(plot_file, transparent=True)
        else:
            plt.show()
        plt.close(fig)


if __name__ == '__main__':
    main()
