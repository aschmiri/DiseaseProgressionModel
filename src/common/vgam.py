#! /usr/bin/env python
# print __doc__
import os.path
import csv
import re
import numpy as np
from src.common import adni_tools as adni

MIN_PROGRESS = -42
MAX_PROGRESS = 51


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
            if not rid in measurements:
                measurements.update({rid: {}})

            # Get scan time
            viscode = row['VISCODE']
            if viscode in measurements[rid]:
                print 'ERROR: Entry already exists:', rid, viscode
                break
            measurements[rid].update({viscode: {}})

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
                print 'ERROR: Invalid diagnosis:', viscode
                break
            measurements[rid][viscode].update({'DX.scan': dx})

            # Get and normalise volumes
            for biomarker in adni.biomarker_names:
                value = adni.safe_cast(row[biomarker])
                # Scale only if biomarker is a volume
                if biomarker in adni.volume_names:
                    value = value / float(row['FactorMNI'])
                measurements[rid][viscode].update({biomarker: value})

    #
    # Add time relative to point of conversion to each data set
    valid_rids = []
    for rid, rid_data in measurements.items():
        data_viscode = []
        data_scantime = []
        data_diagnosis = []

        for viscode, scan_data in rid_data.items():
            if viscode == 'bl':
                scantime = 0
            elif re.match('m[0-9][0-9]', viscode):
                scantime = int(viscode[1:])
            else:
                print 'ERROR: Invalid viscode:', viscode
                break

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
        print 'WARNING: failed to read', pdf_file
        print e
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
        print 'Reading pdfs for', biomarker

        pdf_file = os.path.join(folder, biomarker.replace(' ', '_') + '_densities.csv')
        metric_grid, progress_grid, function_values = _read_pdf_file(pdf_file)

        if metric_grid != None and progress_grid != None and function_values != None:
            pdfs.update({biomarker: {}})
            pdfs[biomarker].update({'values': metric_grid})

            for i in range(len(function_values)):
                pdfs[biomarker].update({progress_grid[i]: function_values[i]})

    return pdfs
