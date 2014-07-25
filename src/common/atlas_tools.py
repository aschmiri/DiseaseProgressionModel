#! /usr/bin/env python
# print __doc__
import numpy as np


################################################################################
#
# read_datafile
#
################################################################################
def find_file(filename, method_folder, type_folder):
    import os.path
    import adni_tools as adni

    folder_ADNI1 = os.path.join(adni.data_folder, 'ADNI1', method_folder, type_folder)
    folder_ADNI2 = os.path.join(adni.data_folder, 'ADNI2', method_folder, type_folder)
    folder_ADNIGO = os.path.join(adni.data_folder, 'ADNIGO', method_folder, type_folder)

    basename = os.path.basename(filename)
    if 'dof' in type_folder:
        basename = basename.replace('.nii.gz', '.dof.gz')
    elif 'images' in type_folder:
        basename = basename.replace('.dof.gz', '.nii.gz')

    filename = os.path.join(folder_ADNI1, basename)

    if os.path.exists(filename):
        return filename

    filename = os.path.join(folder_ADNI2, basename)
    if os.path.exists(filename):
        return filename

    filename = os.path.join(folder_ADNIGO, basename)
    if os.path.exists(filename):
        return filename


################################################################################
#
# read_datafile
#
################################################################################
def read_datafile(datafile, diagnosis='ALL', age_regression=False):
    import csv
    import os.path

    rids = []
    ages = []
    mmses = []
    states = []
    images = []

    with open(datafile, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        headers = reader.next()
        for row in reader:
            rid = int(row[headers.index('RID')])
            dx = row[headers.index('DX.bl')]
            age = float(row[headers.index('AGE')])
            mmse = int(row[headers.index('MMSE')])
            state = float(row[headers.index('DPI')])
            image = row[headers.index('FILE')]
            if os.path.exists(image):
                if diagnosis == 'ALL' or diagnosis in dx or dx in diagnosis:
                    rids.append(rid)
                    ages.append(age)
                    mmses.append(mmse)
                    states.append(state)
                    images.append(image)
    # Cast to np arrays
    states = np.array(states)
    ages = np.array(ages)

    # Perform age regression
    if age_regression:
        import scipy.stats as stats

        mean_state = np.mean(states)
        slope, intercept, _, _, _ = stats.linregress(ages, states)
        states = states - (ages * slope + intercept) + mean_state

    # Return read data
    return np.array(rids), ages, np.array(mmses), states, np.array(images)


################################################################################
#
# adaptive_kernel_regression
#
################################################################################
def adaptive_kernel_regression(states, state,
                               sigma_min=0.05, sigma_max=3.0, sigma_delta=0.05,
                               min_weight=0.01, required_subjects=100):

    for sigma in np.arange(sigma_min, sigma_max + sigma_delta, sigma_delta):
        weights = np.exp(-0.5 * np.square((states - state) / sigma))
        indices = np.where(weights > min_weight)
        if len(indices[0]) >= required_subjects:
            break

    print 'Ended at sigma=' + str(sigma) + ' for state=' + str(state) + ' with ' + str(len(indices[0])) + ' images.'
    return sigma, weights, indices


################################################################################
#
# main
#
################################################################################
if __name__ == "__main__":
    print ' No main implemented.'
