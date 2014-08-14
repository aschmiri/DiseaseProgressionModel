#! /usr/bin/env python2.7
import os.path
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from bmia.common import adni_tools as adni
from bmia.common import adni_plot as aplt


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

    # data_file = os.path.join(adni.project_folder, 'lists/volumes_probabilistic.csv')
    # data_file = os.path.join(adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv')
    # data_file = os.path.join(adni.project_folder, 'lists/volumes_segbased_longitudinal.csv')
    data_file = os.path.join(adni.project_folder, 'lists/volumes_segbased_graphcut.csv')

    volume_name = adni.volume_names[19]
    print 'Analysing', volume_name

    traj_x = []
    traj_y = []
    traj_d = []
    previous_rid = None
    with open(data_file, 'rb') as csvfile:
        rows = csv.DictReader(csvfile)
        for row in rows:
            # Get rid and diagnosis
            rid = int(row['RID'])

            # Get scan time
            viscode = row['VISCODE']
            if viscode == 'bl':
                scan_time = 0
            elif re.match('m[0-9][0-9]', viscode):
                scan_time = int(viscode[1:])
            else:
                print 'ERROR: Invalid viscode:', viscode
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
                print 'ERROR: Invalid diagnosis:', viscode
                break

            # Get and normalise volumes
            volume = float(row[volume_name])
            # volumes = [float(vol) * factor for vol in row]

            if rid != previous_rid and previous_rid is not None:
                # Plot previous subject
                print 'Plotting subject', previous_rid
                plot_trajectories(traj_x, traj_y, traj_d, previous_rid,
                                  biomarker=volume_name)
                traj_x = []
                traj_y = []
                traj_d = []
            traj_x.append(scan_time)
            traj_y.append(volume)
            traj_d.append(dx)

            # Store previous rid for next row
            previous_rid = rid

    print np.mean(list_r)
    print np.mean(list_t0)
    print np.mean(list_up)

    x = np.linspace(0, plot_max, plot_steps)
    y = sigmoid(x, np.mean(list_t0), np.exp(np.mean(list_r)), 1, np.mean(list_up))
    # y = sigmoid(x, np.mean(list_t0), np.mean(list_r), 1, np.mean(list_up))

    plt.title(volume_name + ' mean model ')
    plt.xlabel('Months after baseline')
    plt.ylabel('Volume growth')
    plt.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.5, 0.0)),
                mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.8, 0.8, 0.0)),
                mpl.patches.Rectangle((0, 0), 1, 1, fc=(1.0, 0.0, 0.0))],
               ['CN', 'MCI', 'AD'],
               bbox_to_anchor=(0.25, 0.95))
    plt.plot(x, y, color='b', label='fit1', linewidth=2)
    plt.show()

    for r, up, dx in zip(list_r, list_up, list_dx):
        y = sigmoid(x, plot_max / 2, r, 1, up)
        # z = interpolate_color_scale(xdata, dx, plot_steps)
        # colorline(x, y, z)
        # plt.plot(x, y, color=color)

    plt.show()


def sigmoid(t, t0, r, lower, upper):
    y = lower + (upper - lower) * (1 / (1 + np.exp(-r * (t - t0))))
    return y


def exponential(t, t0, r, r2, lower):
    y = lower + r * np.exp(r2 * t - t0)
    return y


def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    return segments


# Interface to LineCollection:

def colorline(x, y, z=None, cmap=aplt.progression_cmap, norm=plt.Normalize(0.0, 1.0), linewidth=2, alpha=1.0):
    ''' Plot a colored line with coordinates x and y

    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mpl.collections.LineCollection(segments, array=z, cmap=aplt.progression_cmap, norm=norm, linewidth=linewidth, alpha=alpha)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


def interpolate_color_scale(xdata, dx, plot_max, plot_steps):
    x = np.linspace(0, plot_max, plot_steps)

    # Interpolate color values
    ind = 0
    prev_ind = xdata[ind]
    next_ind = xdata[ind + 1]
    prev_val = dx[ind]
    next_val = dx[ind + 1]
    z = np.zeros(plot_steps)
    for i in range(plot_steps):
        if x[i] > next_ind:
            ind += 1
        if ind >= len(dx) - 1:
            z[i] = next_val
        else:
            prev_val = dx[ind]
            next_val = dx[ind + 1]
            prev_ind = xdata[ind]
            next_ind = xdata[ind + 1]
            factor = (x[i] - prev_ind) / (next_ind - prev_ind)
            z[i] = prev_val + (next_val - prev_val) * factor
    return z


def plot_trajectories(traj_x, traj_y, traj_d, rid=0, plot_steps=100, biomarker=''):
    if len(traj_x) > 5:  # and traj_d[-1] == 'CN':

        # Convert and sort data
        xdata = np.array(traj_x)
        ydata = np.array(traj_y)
        ddata = np.array(traj_d)

        args = np.argsort(xdata)
        xdata = xdata[args]
        ydata = ydata[args]
        ddata = ddata[args]

        try:
            # Fit curve
            popt, _ = curve_fit(sigmoid, xdata, ydata, p0=(25, 1, min(ydata), max(ydata)))
            # popt, _ = curve_fit(exponential, xdata, ydata, p0 = (5, 0, 0.01, np.mean(ydata)))
            ydata = ydata

            # if popt[0] > 0 and popt[0] < 100 and popt[3] / popt[2] < 2:
            print 'Fitting result:', popt

            # Plot curve
            x = np.linspace(0, plot_max, plot_steps)
            y = sigmoid(x, *popt)
            z = interpolate_color_scale(xdata, ddata, plot_max, plot_steps)

            colorline(x, y, z)

            list_r.append(np.sign(popt[1]) * np.log(np.abs(popt[1])))
            list_t0.append(popt[0])
            list_up.append(popt[3])
            list_dx.append(ddata[-1])

        except RuntimeError:
            print 'Optimal parameters not found'

        # Plot data
        plt.title(biomarker + ' of subject ' + str(rid))
        plt.xlabel('Months after baseline')
        plt.ylabel('Volume ')
        # plt.axis([-5, 65,  18500, 25000])
        plt.legend([mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.0, 0.5, 0.0)),
                    mpl.patches.Rectangle((0, 0), 1, 1, fc=(0.8, 0.8, 0.0)),
                    mpl.patches.Rectangle((0, 0), 1, 1, fc=(1.0, 0.0, 0.0))],
                   ['CN', 'MCI', 'AD'],
                   bbox_to_anchor=(0.25, 0.95))
        ax = plt.gca()
        ax.scatter(xdata, ydata, c=ddata, cmap=aplt.progression_cmap, vmin=0.0, vmax=1.0, s=40, linewidths=0)
        plt.show()


if __name__ == '__main__':
    main()
