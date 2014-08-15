#! /usr/bin/env python2.7
import os.path
import csv
import math
# import re
import datetime
import numpy as np
from bmia.common import log as log
from bmia.common import adni_tools as adni

MIN_PROGRESS = -1200  # -42
MAX_PROGRESS = 1200  # 51

MIN_DPR = 0.1
MAX_DPR = 3.0


################################################################################
#
# get_measurements_as_collection()
#
################################################################################
def get_measurements_as_collection(data_file):
    '''Return all measurements as a collection.

    Arguments:
    data_file -- the names of the csv file containing the data

    Returns:
    A collection with the following structure:
    { RID : { viscode : { DX.scan : <diagnosis> }
                        { scantime : <months after bl> }
                        { progress : <months relative to conversion> }
                        { <biomarker> : <value> } }
            { viscode : ... }
      RID : ... }
    '''
    measurements = {}

    #
    # Get all measurements from CSV file
    with open(data_file, 'rb') as csvfile:
        rows = csv.DictReader(csvfile)
        for row in rows:
            # Get rid
            rid = int(row['RID'])
            if rid not in measurements:
                measurements.update({rid: {}})

            # Get scan time
            viscode = row['VISCODE']
            if viscode in measurements[rid]:
                print log.ERROR, 'Entry already exists:', rid, viscode
                break
            measurements[rid].update({viscode: {}})

            # Get scan date
            scandate = datetime.datetime.strptime(row['ScanDate'], "%Y-%m-%d").date()
            measurements[rid][viscode].update({'scandate': scandate})

            # Get age
            measurements[rid][viscode].update({'AGE.scan': float(row['AGE.scan'])})

            # Get diagnosis as numerical value
            dx_str = row['DX.scan']
            if dx_str == 'AD':
                dx = 1.0
            elif dx_str == 'MCI':
                dx = 0.5
            elif dx_str == 'CN':
                dx = 0.0
            else:
                print log.ERROR, 'Invalid diagnosis:', viscode
                break
            measurements[rid][viscode].update({'DX.scan': dx})

            # Get and normalise volumes
            for biomarker in adni.biomarker_names:
                if biomarker in row:
                    value = adni.safe_cast(row[biomarker])
                    # Scale only if biomarker is a volume
                    if value > 0 and biomarker in adni.volume_names:
                        value = value / float(row['FactorMNI'])
                    measurements[rid][viscode].update({biomarker: value})

    #
    # Add time relative to point of conversion to each data set
    valid_rids = []
    for rid, rid_data in measurements.items():
        data_viscode = []
        data_scantime = []
        data_diagnosis = []

        if 'bl' not in measurements[rid]:
            print log.WARNING, 'No bl scan for subject {0}!'.format(rid)
            continue

        bl_date = measurements[rid]['bl']['scandate']
        for viscode, scan_data in rid_data.items():
            #             if viscode == 'bl':
            #                 scantime = 0
            #             elif re.match('m[0-9][0-9]', viscode):
            #                 scantime = int(viscode[1:])
            #             else:
            #                 print 'ERROR: Invalid viscode:', viscode
            #                 break
            fu_date = measurements[rid][viscode]['scandate']
            scantime = (fu_date - bl_date).days

            data_viscode.append(viscode)
            data_scantime.append(scantime)
            data_diagnosis.append(scan_data['DX.scan'])

        data_viscode = np.array(data_viscode)
        data_scantime = np.array(data_scantime)
        data_diagnosis = np.array(data_diagnosis)

        args = np.argsort(data_scantime)
        data_viscode = data_viscode[args]
        data_scantime = data_scantime[args]
        data_diagnosis = data_diagnosis[args]

        if data_diagnosis[-1] == 1.0 and data_diagnosis[0] == 0.5:
            valid_rids.append(rid)
            scantime_prev = data_scantime[0]
            for diagnosis, scantime in zip(data_diagnosis, data_scantime):
                if diagnosis == 1.0:
                    time_convert = scantime_prev + (scantime - scantime_prev) / 2
                    break
                else:
                    scantime_prev = scantime

            for viscode, scantime in zip(data_viscode, data_scantime):
                measurements[rid][viscode].update({'scantime': scantime})
                measurements[rid][viscode].update({'progress': scantime - time_convert})

    #
    # Remove non-valid rids
    measurements = {rid: value for rid, value in measurements.items() if rid in valid_rids}

    return measurements


################################################################################
#
# get_rcd_as_collection()
#
################################################################################
def get_rcd_as_collection(measurements):
    '''Return all a collection that indicates for each RID if the subject
    is classified as RCD (rapid cognitive decline.

    Arguments:
    measurements -- the measurements as a collection

    Returns:
    A collection with the following structure:
    { RID : [True|False]
      ... }
    '''
    rcds = {}

    for rid in measurements:
        mmse_bl = measurements[rid]['bl']['MMSE']
        mmse_24 = measurements[rid]['m24']['MMSE']
        # rcd = True if (mmse_24 - mmse_bl) < -7 else False
        decline = mmse_bl - mmse_24
        rcds.update({rid: decline})

    return rcds


################################################################################
#
# _read_pdf_file()
#
################################################################################
def _read_pdf_file(pdf_file):
    function_values = []
    progress_grid = []
    try:
        with open(pdf_file, 'rb') as csvfile:
            rows = csv.reader(csvfile)
            rows.next()
            for row in rows:
                name = row[0]
                data = [float(row[i]) for i in range(1, len(row))]
                if name == 'values':
                    metric_grid = data
                else:
                    progress_grid.append(int(name))
                    function_values.append(data)
        return np.asarray(metric_grid), np.asarray(progress_grid), np.asarray(function_values)
    except Exception as e:
        print log.WARNING, 'Failed to read', pdf_file
        print log.WARNING, e
        return None, None, None


################################################################################
#
# get_pdfs_as_collection()
#
################################################################################
def get_pfds_as_collection(folder=os.path.join(adni.project_folder, 'data'),
                           biomarkers=adni.biomarker_names):
    '''Return all density distribution functions (PDFs) as a collection.

    Keyword arguments:
    biomarkers -- the names of the biomarkers that are to be included
    folder -- the folder where the csv files describing the biomarkers are found

    Returns:
    A collection with the following structure:
    { <biomarker> : { values : [sample points of value] }
                    { MIN_PROGRESS : [PDF at progress MIN_PROGRESS] }
                    ...
                    { MAX_PROGRESS : [PDF at progress MAX_PROGRESS] }
      <biomarker> : ... }
    '''
    pdfs = {}
    for biomarker in biomarkers:
        print log.INFO, 'Reading pdfs for', biomarker

        pdf_file = os.path.join(folder, biomarker.replace(' ', '_') + '_densities.csv')
        metric_grid, progress_grid, function_values = _read_pdf_file(pdf_file)

        if metric_grid is not None and progress_grid is not None and function_values is not None:
            pdfs.update({biomarker: {}})
            pdfs[biomarker].update({'values': metric_grid})

            for i in range(len(function_values)):
                pdfs[biomarker].update({progress_grid[i]: function_values[i]})

    return pdfs


################################################################################
#
# get_scaled_measurements()
#
################################################################################
def get_scaled_measurements(measurements, biomarkers=adni.biomarker_names):
    densities = get_pfds_as_collection(folder=os.path.join(adni.project_folder, 'data'),
                                       biomarkers=biomarkers)
    for rid in measurements:
        print log.INFO, 'Estimating optimal scaling for subject {0}...'.format(rid)

        # Collect samples
        samples = {}
        for viscode in measurements[rid]:
            samples.update({viscode: measurements[rid][viscode]})

        # Get and save scaling
        scaling = get_scaling_for_samples(densities, samples, biomarkers)
        measurements[rid].update({'scaling': scaling})

        # Update all progresses
        for viscode in measurements[rid]:
            if isinstance(viscode, (int, long)):
                progress = measurements[rid][viscode]['progress']
                measurements[rid][viscode].update({'progress': progress * scaling})

    return measurements


################################################################################
#
# get_dpi_for_samples()
#
################################################################################
def get_dpi_for_samples(densities, samples, biomarkers=adni.biomarker_names):
    '''Return the estimated DPI of a subject given a number of samples and
    a set of biomarkers.

    Arguments:
    densities --
    samples --
    biomarkers --

    Returns:
    The estimated DPI
    '''
    max_scantime = np.max([samples[viscode]['scantime'] for viscode in samples])
    dpis = np.arange(MIN_PROGRESS, MAX_PROGRESS - max_scantime, 0.5)
    probs = []
    for dpi in dpis:
        prob = 1.0
        for viscode in samples:
            offset = samples[viscode]['scantime']
            prob *= _get_probability_for_dpi(dpi + offset, densities, samples[viscode], biomarkers=biomarkers)
        probs.append(prob)

    return dpis[np.argmax(probs)]


################################################################################
#
# get_dpi_dpr_for_samples()
#
################################################################################
def get_dpi_dpr_for_samples(densities, samples, biomarkers=adni.biomarker_names):
    '''Return the estimated DPI and DPR of a subject given a number of samples
    and a set of biomarkers.

    Arguments:
    densities --
    samples --
    biomarkers --

    Returns:
    The estimated DPI and DPR
    '''
    max_scantime = np.max([samples[viscode]['scantime'] for viscode in samples])
    dprs = np.arange(MIN_DPR, MAX_DPR, 0.1)
    prob_max = []
    dpi_max = []
    for dpr in dprs:
        dist = max_scantime * dpr
        dpis = np.arange(MIN_PROGRESS, MAX_PROGRESS - dist, 0.5)
        probs = []
        for dpi in dpis:
            prob = 1.0
            for viscode in samples:
                offset = samples[viscode]['scantime'] * dpr
                prob *= _get_probability_for_dpi(dpi + offset, densities, samples[viscode], biomarkers=biomarkers)
            probs.append(prob)

        arg_max = np.argmax(probs)
        prob_max.append(probs[arg_max])
        dpi_max.append(dpis[arg_max])

    # Find the DPR with the highest probability of the corresponding DPI
    arg_max = np.argmax(prob_max)
    return dpi_max[arg_max], dprs[arg_max]


################################################################################
#
# get_scaling_for_samples()
#
################################################################################
def get_scaling_for_samples(densities, samples, biomarkers=adni.biomarker_names):
    '''Return the estimated DPI and DPR of a subject given a number of samples
    and a set of biomarkers.

    Arguments:
    densities --
    samples --
    biomarkers --

    Returns:
    The estimated DPI and DPR
    '''
    scalings = np.arange(0.2, 3, 0.05)
    probs = []
    for scaling in scalings:
        prob = 1.0
        for viscode in samples:
            scaled_progress = samples[viscode]['progress'] * scaling
            prob *= _get_probability_for_dpi(scaled_progress, densities, samples[viscode], biomarkers=biomarkers)
        probs.append(prob)

    # Find the scaling with the highest probability of the corresponding DPI
    arg_max = np.argmax(probs)
    return scalings[arg_max]


################################################################################
#
# _get_probability_for_dpi()
#
################################################################################
def _get_probability_for_dpi(dpi, densities, sample, biomarkers=adni.biomarker_names):
    ''' Get the probability of a DPI given the density functions, a sample and
    a set of biomarkers.

    Arguments:
    dpi --
    densities --
    sample --
    biomarkers --

    Returns
    The probability of the DPI for the given data.
    '''
    prob_sample = 1.0

    prog_prev = math.floor(dpi)
    prog_next = math.ceil(dpi)
    prog_offset = dpi - prog_prev

    for biomarker in biomarkers:
        if biomarker not in densities:
            print log.WARNING, 'No densities available for', biomarker
            prob_sample = 0
        elif prog_prev not in densities[biomarker] or prog_next not in densities[biomarker]:
            # print log.WARNING, 'No densities for time', prog_next
            prob_sample = 0
        else:
            values = densities[biomarker]['values']
            value_sample = sample[biomarker]

            if value_sample is None:
                print log.WARNING, 'Sample has no value for', biomarker
            else:
                # Find value in probability list
                i = 0
                while values[i] < value_sample and i < len(values):
                    i += 1
                i_prev = i - 1 if i > 0 else 0
                i_next = i if i < len(values) else len(values)

                # Get factor for value
                value_prev = values[i_prev]
                value_next = values[i_next]
                factor = (value_sample - value_prev) / (value_next - value_prev)

                # Interpolate probability
                prob_prev = densities[biomarker][prog_prev][i_prev]
                prob_next = densities[biomarker][prog_prev][i_next]
                prob_sample_prev = factor * prob_next + (1 - factor) * prob_prev

                if prog_offset == 0.0:
                    prob_sample *= prob_sample_prev
                else:
                    prob_prev = densities[biomarker][prog_next][i_prev]
                    prob_next = densities[biomarker][prog_next][i_next]
                    prob_sample_next = factor * prob_next + (1 - factor) * prob_prev

                    # Interpolate between
                    prob_sample *= (1 - prog_offset) * prob_sample_prev + prog_offset * prob_sample_next

    return prob_sample


################################################################################
#
# yeojohnson_density()
#
################################################################################
def yeojohnson_density(y, lmbda, mu, sigma):
    '''Return the probability of a value y given lambda, mu and sigma'''
    return (1 / sigma) * \
        _std_normal_dist((_yeojohnson(y, lmbda) - mu) / sigma) * \
        np.power(np.abs(y) + 1, np.sign(y) * (lmbda - 1))


################################################################################
#
# _yeojohnson()
#
################################################################################
def _yeojohnson(y, lmbda):
    if y < 0:
        if lmbda == 2:
            return -np.log(-y + 1)
        else:
            return -(np.power(-y + 1, 2 - lmbda) - 1) / (2 - lmbda)
    else:
        if lmbda == 0:
            return np.log(y + 1)
        else:
            return (np.power(y + 1, lmbda) - 1) / lmbda


################################################################################
#
# _std_normal_dist()
#
################################################################################
def _std_normal_dist(x):
    return np.exp(-0.5 * np.square(x)) / np.sqrt(2 * np.pi)
