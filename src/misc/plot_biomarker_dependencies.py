#! /usr/bin/env python2.7
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import plotting_tools as pt
from common.datahandler import DataHandler


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--biomarkers', nargs=2, default=['D1', 'D2'], help='name of the biomarker to be plotted')
    parser.add_argument('--plot_file', type=str, default=None, help='filename of the output file')
    args = parser.parse_args()

    # Collect data for test
    data_handler = DataHandler.get_data_handler(biomarkers=args.biomarkers)
    biomarkers = data_handler.get_biomarker_names()
    measurements = data_handler.get_measurements_as_dict(biomarkers=biomarkers,
                                                         select_complete=True)

    # Collect biomarker values
    biomarkers_1 = []
    biomarkers_2 = []
    diagnoses = []
    for rid in measurements:
        for visit in measurements[rid]:
            biomarkers_1.append(measurements[rid][visit][biomarkers[0]])
            biomarkers_2.append(measurements[rid][visit][biomarkers[1]])
            diagnoses.append(measurements[rid][visit]['DX.scan'])
    diagnoses = np.array(diagnoses)
    diagnoses[(0.25 <= diagnoses) & (diagnoses <= 0.75)] = 0.5

    # Setup plot
    fig, ax = plt.subplots()
    pt.setup_axes(plt, ax)
    ax.scatter(biomarkers_1, biomarkers_2, s=15.0, c=diagnoses, edgecolor='none',
               vmin=0.0, vmax=1.0, cmap=pt.progression_cmap, alpha=0.25)
    ax.set_xlabel(biomarkers[0])
    ax.set_ylabel(biomarkers[1])

    # Plot legend
    rects = [mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_cn + (0.25,), linewidth=0),
             mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_mci + (0.25,), linewidth=0),
             mpl.patches.Rectangle((0, 0), 1, 1, fc=pt.color_ad + (0.25,), linewidth=0)]
    labels = ['CN', 'MCI', 'AD']
    legend = ax.legend(rects, labels, fontsize=10, ncol=len(rects), loc='upper center', framealpha=0.9)
    legend.get_frame().set_edgecolor((0.6, 0.6, 0.6))

    # Draw or save the plot
    plt.tight_layout()
    if args.plot_file is not None:
        plt.savefig(args.plot_file, transparent=True)
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
