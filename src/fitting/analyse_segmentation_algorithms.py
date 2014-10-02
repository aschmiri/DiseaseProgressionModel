#! /usr/bin/env python2.7
import os.path
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
from common import log as log
from common import adni_tools as adni


def main():
    global list_r
    global list_t0
    global list_up
    global list_dx
    global plot_max
    global plot_steps

    list_r = []
    list_t0 = []
    list_up = []
    list_dx = []
    plot_steps = 100
    plot_max = 60

    data_file_r = os.path.join(adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv')
    data_file_l = os.path.join(adni.project_folder, 'lists/volumes_segbased_longitudinal.csv')
    data_file_g = os.path.join(adni.project_folder, 'lists/volumes_segbased_graphcut.csv')
    data_file_c = os.path.join(adni.project_folder, 'lists/volumes_segbased_consistent.csv')

#     for vol_index in [18, 19, 22, 23]:
    for structure in adni.structure_names:
        vol_index = adni.structure_names.index(structure)

        print log.INFO, 'Plotting {0}...'.format(adni.structure_names[vol_index])

        for rid in [112, 307, 408]:
            plot_points_for_subject(data_file_r, rid, vol_index, 'reg-based', 'r')
            plot_points_for_subject(data_file_l, rid, vol_index, 'longitudinal', 'g')
            plot_points_for_subject(data_file_c, rid, vol_index, 'consistent', 'k')
            plot_points_for_subject(data_file_g, rid, vol_index, 'graphcut', 'b')

            plt.title('{0} of subject {1}'.format(adni.structure_names[vol_index], rid))
            plt.xlabel('Months after baseline')
            plt.ylabel('Volume')
            plt.legend()
            plt.show()


def plot_points_for_subject(data_file, rid, vol_index, label, color):
    print log.INFO, 'Plotting subject {0}...'.format(rid)
    traj_x = []
    traj_y = []
    traj_d = []
    with open(data_file, 'rb') as csvfile:
        rows = csv.DictReader(csvfile)
        for row in rows:
            # Get rid and diagnosis
            if rid == int(row['RID']):
                # Get scan time
                viscode = row['VISCODE']
                if viscode == 'bl':
                    scan_time = 0
                elif re.match('m[0-9][0-9]', viscode):
                    scan_time = int(viscode[1:])
                else:
                    print log.ERROR, 'Invalid viscode: {0}'.format(viscode)
                    break

                # Get diagnosis
                dx_str = row['DX.scan']
                if dx_str == 'AD':
                    dx = 1.0
                elif dx_str == 'MCI':
                    dx = 0.5
                elif dx_str == 'CN':
                    dx = 0.0
                else:
                    print log.ERROR, 'Invalid diagnosis: {0}'.format(dx_str)
                    break

                # Get and normalise volumes
                volume = float(row[adni.structure_names[vol_index]])
                # volumes = [float(vol) * factor for vol in row]

                traj_x.append(scan_time)
                traj_y.append(volume)
                traj_d.append(dx)

    if len(traj_x) > 0:  # and traj_d[-1] == 'CN':

        # Convert and sort data
        xdata = np.array(traj_x)
        ydata = np.array(traj_y)
        ddata = np.array(traj_d)

        args = np.argsort(xdata)
        xdata = xdata[args]
        ydata = ydata[args]
        ddata = ddata[args]

        # Plot data
        plt.plot(xdata, ydata, c=color, label=label)


if __name__ == '__main__':
    main()
