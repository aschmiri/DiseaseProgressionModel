#!/usr/bin/env python
import os
from os.path import join, normpath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import csv
import math
import multiprocessing


def erf(x):
    """ Calculate the error function of the given value. """
    # constants
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    # Save the sign of x
    sign = 1
    if x < 0:
        sign = -1
        x = abs(x)

    # A & S 7.1.26
    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)

    return sign * y


def cdf(z):
    return .5 * (1.0 + erf(z / math.sqrt(2)))


def psi(y, L):
    if y > 0 and L != 0:
        return (((y + 1) ** L) - 1) / L
    elif y > 0 and L == 0:
        return math.log(y + 1)
    elif y < 0 and L != 2:
        return -(((-y + 1) ** (2 - L)) - 1) / (2 - L)
    elif y < 0 and L == 2:
        return -math.log(-y + 1)
    else:
        return None


def z_score(L, M, S, vol, vol_offset):
    return (psi(vol + vol_offset, L) - M) / S


def _parse_output_filename(filename):
    name = os.path.basename(filename).strip('.csv').split('_')
    return {'ref': name[0], 'correction': name[1], 'gender': name[2], 'region': name[3], 'hemisphere': name[4], 'tissue': name[5]}


def _get_reference_filename(curve_path, d):
    filename = os.path.join(curve_path, "_".join([d['ref'], d['correction'], d['gender'], d['region'], d['hemisphere'], d['tissue']]) + '_curves.csv')
    #print "filename: %s" % filename
    if os.path.isfile(filename):
        return filename
    else:
        return None


def _get_z_score(curve_file, patient_age, patient_volume):
    data = csv.reader(open(curve_file, 'rb'), dialect='excel')
    data.next()
    # find the L,M,S and offset for the closest age
    closest = min((r for r in data), key=lambda x: abs(float(x[2]) - patient_age))
    [L, M, S, vol_offset] = map(float, closest[3:6] + [closest[11]])
    return z_score(L, M, S, patient_volume, vol_offset)


def _get_p_value(curve_file, patient_age, patient_volume):
    return 100 * cdf(_get_z_score(curve_file, patient_age, patient_volume))


def plot(output_dir, curve_file, meas_file, plot_points):
    # Load csv in a structured array
    r = mlab.csv2rec(curve_file)
    n = _parse_output_filename(curve_file)
    output_filename = normpath(join(output_dir, ("_".join(["%s"] * 6) % (n['ref'], n['correction'], n['gender'], n['region'], n['hemisphere'], n['tissue'])) + ".png"))

    headers = csv.reader(open(curve_file, 'rb'), dialect='excel').next()

    if plot_points:
        # Collect all measurement points from the measurement file
        m = mlab.csv2rec(meas_file)
        if n['gender'] == 'male':
            m = np.delete(m, np.where(m['gender'] == 'female'), 0)
            pointcolor='black'
        else:
            m = np.delete(m, np.where(m['gender'] == 'male'), 0)
            pointcolor='magenta'
        age_points = m['age']
        vol_points = m['_'.join([n['region'], n['hemisphere'], n['tissue']])]

    # Sort the array on age
    r_sorted = np.sort(r, order=r.dtype.names[2])

    age = r_sorted[r.dtype.names[2]]
    x5 = r_sorted['x5']
    x25 = r_sorted['x25']
    x50 = r_sorted['x50']
    x75 = r_sorted['x75']
    x95 = r_sorted['x95']

    #p_value = _get_p_value(curve_file,patient_age,patient_volume)
    #age_range = age.max() - age.min()
    #print "age range: %0.2f" % age_range

    linestyle = ['g-', 'g-', 'b-', 'g-', 'g-']

    #fig = plt.figure(figsize=(8, 5), dpi=100)
    fig = plt.figure(figsize=(13, 7), dpi=100)

    # Plot the percintile curves
    for (curve, linestyle) in zip([x5, x25, x50, x75, x95], linestyle):
        plt.plot(list(age), list(curve), linestyle, linewidth=2.0)
    if plot_points:
        point, = plt.plot(age_points, vol_points, 'mo', mfc='none', mec=pointcolor, markersize=4)

    x_offset = 0.3
    y_offset = 0.1
    fs = 20
    fw = 'bold'
    # Plot text next to the percentile curves
    for (curve, label) in zip([x5, x25, x50, x75, x95], ['5%', '25%', '50%', '75%', '95%']):
        plt.text(age[-1] + x_offset, curve[-1] - y_offset, label, fontsize=18)
    #plt.text( patient_age+2*x_offset,patient_volume-y_offset, ('age=%0.1f\nvol=%0.1f\np=%i%%' % (patient_age,patient_volume,p_value)), bbox={'facecolor':'#cccccc','edgecolor':'#cccccc','alpha':0.5,'pad':5} )
    #plt.text( patient_age+3*x_offset,patient_volume-y_offset, ('p=%i%%' % (p_value)), bbox={'facecolor':'#cccccc','edgecolor':'#cccccc','alpha':0.5,'pad':5}, fontsize=fs )

    plt.title('Percentile curves for ' + (" ".join(["%s"] * 5) % (n['correction'], n['gender'], n['region'], n['hemisphere'], n['tissue'])), fontsize=fs, fontweight=fw)
    plt.xlabel('Age [y]', fontsize=fs, fontweight=fw)
    plt.ylabel('Volume %s [%s%s]' % (n['correction'], n['tissue'] == 'wml' and 'ln ' or '', n['correction'] == 'abs' and 'ml' or '%'), fontsize=fs, fontweight=fw)
    plt.grid()
    plt.setp(plt.gca().get_xticklabels(), fontsize=18)
    plt.setp(plt.gca().get_yticklabels(), fontsize=18)
    
    if n['region'] == 'all' and n['tissue'] == 'gm':  
        print 'y-scale is set manually.'
        #plt.set_autoscaley_on(False)
        plt.ylim([44,52]) #all total GM
    if n['region'] == 'all' and n['tissue'] == 'wm':
        print 'y-scale is set manually.'
        #plt.set_autoscaley_on(False)
        plt.ylim([28,40]) #all total WM
    if n['region'] == 'fl' and n['tissue'] == 'wm':
        print 'y-scale is set manually.'
        #plt.set_autoscaley_on(False)
        plt.ylim([9,16]) #fl total WM
    if n['region'] == 'fl' and n['tissue'] == 'wml': 
        print 'y-scale is set manually.'
        #plt.set_autoscaley_on(False)
        plt.ylim([-6,2]) #fl total WML

    plt.tight_layout()

    print output_filename
    # Draw the plot
    plt.savefig(output_filename, dpi=100)
    #plt.show()


def chunks(l, n):
    """ Yield succesive n-sized chunks from l. """
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def plot_worker(args):
    output_dir, curve_file, meas_file, plot_points = args
    plot(output_dir, curve_file, meas_file, plot_points)


if __name__ == "__main__":
    # a point is a list of tuples: [(patient_age,patient_volume)]
    #curve_file = '/media/mkoek/DATAPART1/2010_Care4Me/ref/RefDatCurves/ref_ptbv_male_all_total_gm_curves.csv'

    meas_file =  '/media/Data/source/repositoriesBitBucket/bigr-reference-curves/ref_picv.csv'
    curve_path = '/media/Data/source/repositoriesBitBucket/bigr-reference-curves/RefDatRCurves'
    output_dir = '/media/Data/source/repositoriesBitBucket/bigr-reference-curves/output'

    plot_points = False

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    plot_args = []
    for root, path, files in os.walk(curve_path):
        for f in [f for f in files if f.endswith('.csv')]:
            curve_file = join(root, f)
            plot_args += [(output_dir, curve_file, meas_file, plot_points)]
    pool_size = multiprocessing.cpu_count() * 2
    pool = multiprocessing.Pool(processes=pool_size)
    result = pool.map(plot_worker, plot_args)
