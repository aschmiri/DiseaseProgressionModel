#! /usr/bin/env python2.7
import argparse
import numpy as np
import scipy.stats as stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from vgam.datahandler import DataHandler
import fitting.vgam_evaluation as ve


def main():
    parser = argparse.ArgumentParser()
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('-p', '--predict_biomarker', type=str, default='MMSE', help='the biomarker to predict')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output file')
    args = parser.parse_args()

    visits = ['bl', 'm12', 'm24']
    methods = ['cog', 'long', 'mbl', 'img', 'all']
    values = {}
    for method in methods:
        values.update({method: {}})
        _, values_observed, values_naive, values_model = \
            ve.get_biomarker_predictions(args, visits, args.predict_biomarker,
                                         method=method, estimate_dprs=False,
                                         exclude_cn=True, consistent_data=True)
        values[method].update({'observed': values_observed})
        values[method].update({'naive': values_naive})
        values[method].update({'model_dpi': values_model})

        _, values_observed, values_naive, values_model = \
            ve.get_biomarker_predictions(args, visits, args.predict_biomarker,
                                         method=method, estimate_dprs=True,
                                         exclude_cn=True, consistent_data=True)
        values[method].update({'model_dpi_dpr': values_model})

    plot_errors(args, values, methods)


def plot_errors(args, values, methods):
    fig, ax = plt.subplots()
    ve.setup_axes(plt, ax, xgrid=False)

    num_methods = len(methods)
    for i in xrange(num_methods):
        ax.axvline(i * 2 + 1.5, linestyle='-', color='k', alpha=0.4)

    ax.set_axisbelow(True)
    ax.set_title('Prediction of {0}'.format(args.predict_biomarker))
    ax.set_ylabel('Prediction error')
    ax.set_xticklabels([])

    # Append plot for naive estimation
    values_observed = np.array(values['all']['observed'])
    values_naive = np.array(values['all']['naive'])
    errors_naive = np.abs(values_observed - values_naive)
    data = [errors_naive]
    medians = [np.median(errors_naive)]
    weights = ['normal']
    for method in methods:
        values_model1 = np.array(values[method]['model_dpi'])
        values_model2 = np.array(values[method]['model_dpi_dpr'])
        errors_model1 = np.abs(values_observed - values_model1)
        errors_model2 = np.abs(values_observed - values_model2)

        data.append(errors_model1)
        data.append(errors_model2)

        medians.append(np.median(errors_model1))
        medians.append(np.median(errors_model2))

        weights.append('bold' if stats.ttest_rel(errors_naive, errors_model1)[1] < 0.01 else 'normal')
        weights.append('bold' if stats.ttest_rel(errors_naive, errors_model2)[1] < 0.01 else 'normal')

    # Set colors
    colors = [(0.0, 0.0, 0.0)]
    for i in range(num_methods):
        colors.append((0.0, 0.0, 0.5))
        colors.append((0.0, 0.5, 0.0))

    # Set limits
    max_data = np.max(data)
    ax.set_ylim(-0.01 * max_data, 1.01 * max_data)

    # Add labels to x axis
    for i, method in enumerate(methods):
        ax.text(i * 2 + 2.5, -0.04 * max_data, method, horizontalalignment='center')

    # Write median errors as text
    pos = np.arange(1, 2 * num_methods + 2)
    upper_labels = [str(np.round(m, 2)) for m in medians]
    for i in range(2 * num_methods + 1):
        ax.text(pos[i], 0.91 * max_data, upper_labels[i], weight=weights[i],
                horizontalalignment='center', size=10, color=colors[i])

    # Plot boxplots
    boxplot = plt.boxplot(data, patch_artist=True)
    for i in range(len(boxplot['boxes'])):
        ve.set_boxplot_color(boxplot, i, colors[i])

    # Plot legend
    legend = ax.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.0, 0.2), ec=(0.0, 0.0, 0.0, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.0, 0.5, 0.2), ec=(0.0, 0.0, 0.5, 1.0), linewidth=1),
                        mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.5, 0.0, 0.2), ec=(0.0, 0.5, 0.0, 1.0), linewidth=1)],
                       ['Naive', 'DPI', 'DPI + DPR'], fontsize=10, ncol=3, loc='upper center', framealpha=0.9)
    legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    # Show or save plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
