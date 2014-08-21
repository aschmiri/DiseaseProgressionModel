#! /usr/bin/env python2.7
import os.path
import argparse
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cmx
from common import log as log
from common import adni_tools as adni
from common import adni_plot as aplt
from common import vgam as vgam


def main():
    parser = argparse.ArgumentParser()
    parser = vgam.add_common_arguments(parser)
    parser.add_argument('-e', '--extrapolator', type=str, choices=['lin', 'sqrt', 'exp'], default='sqrt', help='the type of extrapolator')
    parser.add_argument('-p', '--no_points', action='store_true', default=False, help='do not plot points')
    parser.add_argument('-m', '--no_mu', action='store_true', default=False, help='do not plot mu')
    parser.add_argument('-d', '--no_densities', action='store_true', default=False, help='do not plot densities')
    parser.add_argument('-s', '--save_file', action='store_true', default=False, help='save the plots as a file')
    args = parser.parse_args()

    biomarker_set = vgam.get_biomarker_set(args)
    data_folders = vgam.get_data_folders(args)

    for biomarker in biomarker_set:
        data_folder = vgam.get_data_folder(data_folders, biomarker)
        plot_model(args, biomarker, data_folder)


def plot_model(args, biomarker, data_folder):
    points_file = os.path.join(data_folder, biomarker.replace(' ', '_') + '.csv')
    curves_file = points_file.replace('.csv', '_curves.csv')
    if not os.path.isfile(points_file) or not os.path.isfile(curves_file):
        print log.ERROR, 'Files not available for {0}!'.format(biomarker)
        return

    print log.INFO, 'Generating plot for {0}...'.format(biomarker)

    #
    # Read model
    #
    pm = ProgressionModel(curves_file, extrapolator=args.extrapolator)
    min_progression_extra = pm.min_progression - 0.3 * (pm.max_progression - pm.min_progression)
    max_progression_extra = pm.max_progression + 0.3 * (pm.max_progression - pm.min_progression)
    progression_linspace = np.linspace(min_progression_extra, max_progression_extra, 100)

    #
    # Setup plot
    #
    plt.figure(figsize=(12, 5), dpi=100)
    if not args.no_densities:
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
    else:
        ax1 = plt.subplot(1, 1, 1)

    ax1.set_title('Percentile curves for {0}'.format(biomarker))
    ax1.set_xlabel('Disease progression relative to point of conversion')
    ax1.set_ylabel('Volume' if biomarker in adni.volume_names else 'Score')
    ax1.set_xlim(min_progression_extra, max_progression_extra)
    ax1.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.8, 0.8, 0.0), linewidth=0),
                mpl.patches.Rectangle((0, 0), 1, 1, fc=(1.0, 0.0, 0.0), linewidth=0)],
               ['MCI', 'AD'], fontsize=10)

    #
    # Plot the percentile curves
    #
    ax1.axvline(pm.min_progression, color='0.15', linestyle=':', alpha=0.8)
    ax1.axvline(pm.max_progression, color='0.15', linestyle=':', alpha=0.8)

    stds = [-1.5, -1.0, -0.5, 0, +0.5, +1.0, +1.5]
    greyvals = ['0.45', '0.3', '0.15', '0', '0.15', '0.3', '0.45']
    min_vals = []
    max_vals = []
    for (greyval, std) in zip(greyvals, stds):
        curve = []
        for progr in progression_linspace:
            lmbda, mu, sigma, yoffset = pm.get_parameters(progr)
            curve.append(vgam.yeojohnson_quantile(lmbda, mu, sigma, std) - yoffset)
        min_vals.append(np.min(curve))
        max_vals.append(np.max(curve))
        ax1.plot(progression_linspace, curve, color=greyval)

        label = '${0} \sigma$'.format(std)
        ax1.text(pm.progressions[-1], curve[-1], label, fontsize=11)

    min_val = np.min(min_vals)
    max_val = np.max(max_vals)

    #
    # Plot parameters
    #
    if not args.no_mu:
        # Get second axis of plot 1
        ax1b = ax1.twinx()

        # Plot all progressions
        ax1b.scatter(pm.all_progressions, pm.all_mus, marker='o', linewidths=0, alpha=0.2)
        ax1b.text(pm.progressions[-1], pm.mus[-1], '$\mu$', color='b', fontsize=11)

        # Plot binned progressions
        ax1b.scatter(pm.progressions, pm.mus, color='b', marker='x')

        # Plot interpolated model
        mus = [pm.get_mu(p) for p in progression_linspace]
        ax1b.plot(progression_linspace, mus, color='b')
        ax1b.set_xlim(min_progression_extra, max_progression_extra)

    #
    # Plot points
    #
    if not args.no_points:
        m = mlab.csv2rec(points_file)
        progr_points = m['progress']
        value_points = m['value']
        diagn_points = [0.5 if p < 0 else 1.0 for p in progr_points]

        print log.INFO, 'Plotting {0} sample points...'.format(len(progr_points))
        ax1.scatter(progr_points, value_points, c=diagn_points, vmin=0.0, vmax=1.0, linewidths=0, cmap=aplt.progression_cmap, alpha=0.25)

    #
    # Plot PDFs
    #
    if not args.no_densities:
        ax2.set_title('Probability density function for {0}'.format(biomarker))
        ax2.set_xlabel('Volume' if biomarker in adni.volume_names else 'Score')
        ax2.set_ylabel('Probability')

        values = np.linspace(min_val, max_val, 250)
        progr_samples = [-2000, -1500, -1000, -500, 0, 500, 1000, 1500, 2000]

        sample_cmap = cmx.ScalarMappable(
            norm=colors.Normalize(vmin=-len(progr_samples) + 1, vmax=(len(progr_samples) - 1)),
            cmap=plt.get_cmap(aplt.progression_cmap))

        for progr in progr_samples:
            sample_color = sample_cmap.to_rgba(progr_samples.index(progr))
            ax1.axvline(progr, color=sample_color, linestyle='--', alpha=0.3)
            ax2.set_xlim(min_val, max_val)

            linestyle = '--' if progr < pm.min_progression or progr > pm.max_progression else '-'
            lmbda, mu, sigma, yoffset = pm.get_parameters(progr)
            probs = [vgam.yeojohnson_density(yoffset + v, lmbda, mu, sigma) for v in values]
            ax2.plot(values, probs, label=str(progr), color=sample_color, linestyle=linestyle)

        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(handles2, labels2, fontsize=10)

    #
    # Draw or save the plot
    #
    plt.tight_layout()
    if args.save_file:
        plot_filename = curves_file.replace('.csv', '.pdf')
        plt.savefig(plot_filename, dpi=100)
    else:
        plt.show()


class ProgressionModel:
    def __init__(self, model_file, extrapolator='sqrt'):
        self.extrapolator = extrapolator
        self.__intitialze_model(model_file)

    def __intitialze_model(self, model_file):
        r = mlab.csv2rec(model_file)
        r_sorted = np.sort(r, order=r.dtype.names[2])
        progrs = r_sorted[r.dtype.names[2]]

        parameters = {}
        for i, progr in enumerate(progrs):
            if progr not in parameters:
                parameters.update({progr: {}})
                parameters[progr].update({'sigma': np.exp(r_sorted['logsigma'][i])})
                parameters[progr].update({'lambda': r_sorted['lambda'][i]})
                parameters[progr].update({'mu': r_sorted['mu'][i]})
            parameters[progr].update({'yoffset': r_sorted['fitmiscyoffset'][i]})

        self.all_progressions = sorted(parameters.keys())
        self.sigma = parameters[self.all_progressions[0]]['sigma']
        self.lmbda = parameters[self.all_progressions[0]]['lambda']
        self.yoffset = parameters[self.all_progressions[0]]['yoffset']

        self.all_mus = [parameters[progr]['mu'] for progr in self.all_progressions]

        self.progressions, self.mus = self.__compute_binned_data(self.all_progressions, self.all_mus)

        self.min_progression = np.min(self.progressions)
        self.max_progression = np.max(self.progressions)

    def get_parameters(self, progression):
        return self.lmbda, self.get_mu(progression), self.sigma, self.yoffset

    def get_mu(self, progression):
        if progression in self.progressions:
            return self.mus[self.progressions.index(progression)]

        elif progression < self.progressions[0]:
            delta_p = self.progressions[1] - self.progressions[0]
            delta_m = self.mus[1] - self.mus[0]

            if self.extrapolator == 'lin':
                # Linear:
                return self.mus[1] - delta_m * ((self.progressions[1] - progression) / delta_p)
            elif self.extrapolator == 'sqrt':
                # Square root:
                return self.mus[1] - delta_m * math.sqrt((self.progressions[1] - progression) / delta_p)
            elif self.extrapolator == 'exp':
                # Exponential:
                if delta_m > 0:
                    return self.mus[0] - 1 + math.exp(delta_m / delta_p * (progression - self.progressions[0]))
                else:
                    return self.mus[0] + 1 - math.exp(-delta_m / delta_p * (progression - self.progressions[0]))
            else:
                print log.ERROR, 'Unknown extrapolator!'
                return 0

        elif progression > self.progressions[-1]:
            delta_p = self.progressions[-1] - self.progressions[-2]
            delta_m = self.mus[-1] - self.mus[-2]

            if self.extrapolator == 'lin':
                # Linear:
                return self.mus[-2] + delta_m * ((progression - self.progressions[-2]) / delta_p)
            elif self.extrapolator == 'sqrt':
                # Square root:
                return self.mus[-2] + delta_m * math.sqrt((progression - self.progressions[-2]) / delta_p)
            elif self.extrapolator == 'exp':
                # Exponential:
                if delta_m > 0:
                    return self.mus[-1] + 1 - math.exp(-delta_m / delta_p * (progression - self.progressions[-1]))
                else:
                    return self.mus[-1] - 1 + math.exp(delta_m / delta_p * (progression - self.progressions[-1]))
            else:
                print log.ERROR, 'Unknown extrapolator!'
                return 0

        else:
            idx_next = 1
            while self.progressions[idx_next] < progression:
                    idx_next += 1
            idx_prev = idx_next - 1

            progression_prev = self.progressions[idx_prev]
            progression_next = self.progressions[idx_next]
            factor = (progression - progression_prev) / (progression_next - progression_prev)

            mu_prev = self.mus[idx_prev]
            mu_next = self.mus[idx_next]
            return (1 - factor) * mu_prev + factor * mu_next

    def __compute_binned_data(self, data_x, data_y, bin_size=92):
        bins_x = []
        bins_y = []

        for curr_x in range(np.min(data_x), np.max(data_x), bin_size):
            bin_x = []
            bin_y = []
            for x, y in zip(data_x, data_y):
                if x >= curr_x and x < curr_x + bin_size:
                    bin_x.append(x)
                    bin_y.append(y)

            if len(bin_x) > 1:
                selected_x = np.array(bin_x)
                selected_y = np.array(bin_y)

                bins_x.append(np.mean(selected_x))
                bins_y.append(np.mean(selected_y))

        return bins_x, bins_y


if __name__ == '__main__':
    main()
