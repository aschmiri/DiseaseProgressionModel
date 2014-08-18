#! /usr/bin/env python2.7
import os.path
import csv
import math
import datetime
import numpy as np
from common import log as log
from common import adni_tools as adni

MIN_PROGRESS = -1200  # -42
MAX_PROGRESS = 1200  # 51

MIN_DPR = 0.1
MAX_DPR = 3.0


################################################################################
#
# get_measurements_as_collection()
#
################################################################################
def get_measurements_as_collection(select_converters=True, no_regression=False,
                                   data_file_meta=os.path.join(adni.project_folder, 'lists/metadata.csv'),
                                   data_file_cog=os.path.join(adni.project_folder, 'lists/metadata.csv'),
                                   data_file_vol=os.path.join(adni.project_folder, 'lists/volumes_graphcut.csv'),
                                   data_file_mbl=os.path.join(adni.project_folder, 'lists/manifold_features.csv')):
    '''Return all subjects measurements as a collection.

    Arguments:
    select_converters -- only select MCI -> AD converters
    no_regression -- do not perform age regression
    data_file_meta -- the name of the csv file containing the metadata
    data_file_cog -- the name of the csv file containing the cognitive scores
    data_file_vol -- the name of the csv file containing the structure volumes
    data_file_mbl -- the name of the csv file containing the MBL coordinates

    Returns:
    A collection with the following structure:
    { <rid> : { <viscode> : { DX.scan : <diagnosis> }
                            { AGE.scan : <age in years> }
                            { scandate : <date of scan> }
                            { scantime : <days after bl> }
                            { progress : <days relative to conversion> }
                            { <biomarker1> : <volume> }
                            ... }
              { <viscode> : ... }
      <rid> : ... }
    '''
    # Read data from lists
    metadata = _get_metadata_as_collection(data_file_meta, select_converters=select_converters)

    if data_file_cog is not None:
        metadata = _update_metadata_with_biomarker_values(metadata, data_file_cog,
                                                          adni.cog_score_names,
                                                          no_regression=no_regression)

    if data_file_vol is not None:
        metadata = _update_metadata_with_biomarker_values(metadata, data_file_vol,
                                                          adni.volume_names,
                                                          no_regression=no_regression)

    if data_file_mbl is not None:
        metadata = _update_metadata_with_biomarker_values(metadata, data_file_mbl,
                                                          adni.manifold_coordinate_names,
                                                          no_regression=no_regression)

    # Return measurements
    return metadata


################################################################################
#
# _get_metadata_as_collection()
#
################################################################################
def _get_metadata_as_collection(data_file, select_converters=True):
    '''Return all subjects metadata as a collection.

    Arguments:
    data_file -- the names of the csv file containing the data
    select_converters -- only select MCI -> AD converters

    Returns:
    A collection with the following structure:
    { <rid> : { <viscode> : { DX.scan : <diagnosis> }
                            { AGE.scan : <age in years> }
                            { scandate : <date of scan> }
                            { scantime : <days after bl> }
                            { progress : <days relative to conversion> } }
              { <viscode> : ... }
      <rid> : ... }
    '''
    print log.INFO, 'Collecting metadata...'
    metadata = {}

    if not os.path.isfile(data_file):
        print log.ERROR, 'Data file dies not exist:', data_file
        return metadata

    #
    # Get all measurements from CSV file
    with open(data_file, 'rb') as csvfile:
        rows = csv.DictReader(csvfile)
        for row in rows:
            # Get rid
            rid = int(row['RID'])
            if rid not in metadata:
                metadata.update({rid: {}})

            # Get scan time
            viscode = row['VISCODE']
            if viscode in metadata[rid]:
                print log.WARNING, 'Entry already exists {0} ({1}). Skipping.'.format(rid, viscode)
                break
            metadata[rid].update({viscode: {}})

            # Get scan date
            scandate = datetime.datetime.strptime(row['ScanDate'], "%Y-%m-%d").date()
            metadata[rid][viscode].update({'scandate': scandate})

            # Get age
            metadata[rid][viscode].update({'AGE.scan': float(row['AGE.scan'])})

            # Get factor
            metadata[rid][viscode].update({'FactorMNI': float(row['FactorMNI'])})

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
            metadata[rid][viscode].update({'DX.scan': dx})

    #
    # Add time relative to point of conversion to each data set
    if select_converters:
        print log.INFO, 'Selecting converters...'
        valid_rids = []
        for rid, rid_data in metadata.items():
            data_viscode = []
            data_scantime = []
            data_diagnosis = []

            if 'bl' not in metadata[rid]:
                print log.WARNING, 'No bl scan for subject {0}!'.format(rid)
                continue

            bl_date = metadata[rid]['bl']['scandate']
            for viscode, scan_data in rid_data.items():
                fu_date = metadata[rid][viscode]['scandate']
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
                    metadata[rid][viscode].update({'scantime': scantime})
                    metadata[rid][viscode].update({'progress': scantime - time_convert})

        #
        # Remove non-valid rids
        metadata = {rid: value for rid, value in metadata.items() if rid in valid_rids}

    return metadata


################################################################################
#
# _update_metadata_with_biomarker_values()
#
################################################################################
def _update_metadata_with_biomarker_values(metadata, data_file,
                                           biomarker_names=adni.biomarker_names,
                                           no_regression=False):
    '''
    Update the metadata collection with the biomarker values.

    Arguments:
    metadata -- the metadata
    data_file -- the data file containing the biomarker values
    biomarker_names -- the biomarkers to be read
    no_regression -- do not perform age regression
    '''
    biomarker_values = _get_biomarker_values_as_collection(data_file, metadata,
                                                           biomarker_names=biomarker_names,
                                                           no_regression=no_regression)

    # Update metadata with feature values
    for rid in metadata:
        for viscode in metadata[rid]:
            for biomarker in biomarker_names:
                try:
                    value = biomarker_values[rid][viscode][biomarker]
                    metadata[rid][viscode].update({biomarker: value})
                except:
                    # TODO: Debug message
                    pass

    return metadata


################################################################################
#
# _get_biomarker_values_as_collection()
#
################################################################################
def _get_biomarker_values_as_collection(data_file, metadata,
                                        biomarker_names=adni.biomarker_names,
                                        no_regression=False):
    '''Return all measurements as a collection.

    Arguments:
    data_file -- the names of the csv file containing the data
    biomarker_names -- the biomarkers to be read
    no_regression -- do not perform age regression

    Returns:
    A collection with the following structure:
    { <rid> : { <viscode> : { <biomarker1> : <value> }
                            { <biomarker2> : <value> } ... }
              { <viscode> : ... }
      <rid> : ... }
    '''
    if data_file is None:
        return {}

    NO_REGRESSION_STR = '{0} no regression'
    print log.INFO, 'Collecting biomarker values...'
    values = {}

    #
    # Get all measurements from CSV file
    with open(data_file, 'rb') as csvfile:
        rows = csv.DictReader(csvfile)
        for row in rows:
            # Get rid
            rid = int(row['RID'])
            if rid not in values:
                values.update({rid: {}})

            # Get scan time
            viscode = row['VISCODE']
            if viscode in values[rid]:
                print log.WARNING, 'Entry already exists {0} ({1}). Skipping.'.format(rid, viscode)
                break
            values[rid].update({viscode: {}})

            # Get and normalise volumes
            for biomarker in biomarker_names:
                if biomarker in row:
                    value = adni.safe_cast(row[biomarker])

                    # Normalise volumes with mni transformation
                    if value > 0 and biomarker in adni.volume_names:
                        try:
                            factor = metadata[rid][viscode]['FactorMNI']
                            value /= factor
                        except:
                            break

                    if no_regression:
                        values[rid][viscode].update({biomarker: value})
                    else:
                        values[rid][viscode].update({NO_REGRESSION_STR.format(biomarker): value})

    if not no_regression:
        print log.INFO, 'Performing age regression..'
        for biomarker in biomarker_names:
            rids = []
            viscodes = []
            vals = []
            ages = []
            for rid, visits in values.items():
                for viscode in visits:
                    # Get age
                    try:
                        age = metadata[rid][viscode]['AGE.scan']
                    except:
                        age = None

                    if age is not None and NO_REGRESSION_STR.format(biomarker) in visits[viscode]:
                        value = visits[viscode][NO_REGRESSION_STR.format(biomarker)]
                        if value is not None:
                            ages.append(age)
                            vals.append(value)
                            rids.append(rid)
                            viscodes.append(viscode)

            if len(ages) > 1:
                regressed_vals = _age_regression(ages, vals)

                for rid, viscode, value in zip(rids, viscodes, regressed_vals):
                    values[rid][viscode].update({biomarker: value})

    return values


################################################################################
#
# _age_regression()
#
################################################################################
def _age_regression(timepoints, values):
    '''
    Perform age regression.
    '''
    timepoints = np.array(timepoints)
    values = np.array(values)
    mean_val = np.mean(values)

    import scipy.stats as stats
    slope, intercept, _, _, _ = stats.linregress(timepoints, values)
    values = values - (timepoints * slope + intercept) + mean_val

    return values


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
    { <rid> : [True|False]
      ... }
    '''
    rcds = {}

    for rid in measurements:
        try:
            mmse_bl = measurements[rid]['bl']['MMSE']
            mmse_24 = measurements[rid]['m24']['MMSE']
            # rcd = True if (mmse_24 - mmse_bl) < -7 else False
            decline = mmse_bl - mmse_24
            rcds.update({rid: decline})
        except Exception:
            # TODO: debug message
            pass

    return rcds


################################################################################
#
# get_pdfs_as_collection()
#
################################################################################
def get_pfds_as_collection(folder=os.path.join(adni.project_folder, 'data'),
                           biomarkers=adni.biomarker_names):
    '''Return all density distribution functions (PDFs) as a collection.

    Keyword arguments:
    folder -- the folder where the csv files describing the biomarkers are found
    biomarkers -- the names of the biomarkers that are to be included

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
def interpolate_pdf(pdfs, progression):
    '''
    Interpolate the PDF for a certain progression.
    '''
    if progression in pdfs:
        return pdfs[progression]
    else:
        progressions = pdfs.keys()
        progressions = np.sort(np.array(progressions))

        if progression < progressions[0]:
            print log.WARNING, 'Progression out of scope, returning smallest PDF.'
            return pdfs[progressions[0]]

        elif progression > progressions[-1]:
            print log.WARNING, 'Progression out of scope, returning largest PDF.'
            return pdfs[progressions[-1]]

        else:
            idx = 0
            while progressions[idx] < progression:
                idx += 1

            factor = (progressions[idx] - progression) / (progressions[idx] - progressions[idx - 1])
            pdf_smaller = np.array(pdfs[progressions[idx - 1]])
            pdf_larger = np.array(pdfs[progressions[idx]])

            return factor * pdf_smaller + (1 - factor) * pdf_larger


################################################################################
#
# get_scaled_measurements()
#
################################################################################
def get_scaled_measurements(measurements,
                            folder=os.path.join(adni.project_folder, 'data'),
                            biomarkers=adni.biomarker_names):
    densities = get_pfds_as_collection(folder, biomarkers=biomarkers)
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
    '''Return the optimal scaling value for a subject given a number of samples
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
