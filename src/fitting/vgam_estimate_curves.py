#!/usr/bin/env python

import argparse
import os.path
import csv
import re
import numpy as np
import common.adni_tools as adni
import subprocess
import joblib as jl

def collect_measurements():    
    measurements = {}
    
    #
    # Get all measurements from CSV file
    data_file = os.path.join( adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv' )
    with open( data_file, 'rb' ) as csvfile:
        rows = csv.DictReader( csvfile )
        for row in rows:
            # Get rid
            rid = int( row['RID'] )
            if not rid in measurements:
                measurements.update( {rid : {}} )
            
            # Get scan time
            viscode = row['VISCODE']
            if viscode in measurements[rid]:
                print 'ERROR: Entry already exists:', rid, viscode
                break
            measurements[rid].update( {viscode : {}} )
            
            # Get age
            measurements[rid][viscode].update( {'AGE.scan' : float( row['AGE.scan'] )} )
                         
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
            measurements[rid][viscode].update( {'DX.scan' : dx } )
            
            # Get and normalise volumes
            for score_name in adni.biomarker_names:
                value = adni.safe_cast( row[score_name] )
                if score_name in adni.volume_names:
                    value = value / float( row['FactorMNI'] )
                measurements[rid][viscode].update( {score_name : value } )
    
    #
    # Add time relative to point of conversion to each data set
    valid_rids = []
    for rid, rid_data in measurements.items():
        data_viscode   = []
        data_scantime  = []
        data_diagnosis = []
        
        for viscode, scan_data in rid_data.items():
            if viscode == 'bl':
                scantime = 0
            elif re.match('m[0-9][0-9]', viscode):
                scantime = int( viscode[1:] )
            else:
                print 'ERROR: Invalid viscode:', viscode
                break
            
            data_viscode.append( viscode )
            data_scantime.append( scantime )
            data_diagnosis.append( scan_data['DX.scan'] )
        
        data_viscode   = np.array( data_viscode )
        data_scantime  = np.array( data_scantime )
        data_diagnosis = np.array( data_diagnosis )

        args = np.argsort( data_scantime )
        data_viscode   = data_viscode[args]
        data_scantime  = data_scantime[args]
        data_diagnosis = data_diagnosis[args]
        
        if data_diagnosis[-1] == 1.0 and data_diagnosis[0] == 0.5:
            valid_rids.append( rid )
            scantime_prev = data_scantime[0]
            for diagnosis, scantime in zip(data_diagnosis, data_scantime):
                if diagnosis == 1.0:
                    time_convert = scantime_prev + (scantime-scantime_prev) / 2
                    break
                else:
                    scantime_prev = scantime
    
            for viscode,scantime in zip( data_viscode, data_scantime ):
                measurements[rid][viscode].update( {'scantime' : scantime} )
                measurements[rid][viscode].update( {'progress' : scantime - time_convert} )

    #
    # Remove non-valid rids        
    measurements = { rid : value for rid, value in measurements.items() if rid in valid_rids }
    
    return measurements

def generate_csv_files():
    measurements = collect_measurements()
    for score_name in adni.biomarker_names:
        print 'Generating output CSV for', score_name
        
        csv_file = os.path.join( adni.project_folder, 'data', score_name.replace( ' ', '_' ) + '.csv' )
        writer = csv.writer( open( csv_file, 'wb' ), delimiter=',' )
        writer.writerow(['rid', 'progress', 'value'])

        for rid, rid_data in measurements.items():
            for _, scan_data in rid_data.items():
                progress = adni.safe_cast( scan_data['progress'], int )
                value    = adni.safe_cast( scan_data[score_name], int )
                if progress != None and value != None:
                    writer.writerow( [rid, progress, value] )


def fit_data( args ):
    jl.Parallel( n_jobs=args.nr_threads )( jl.delayed(fit_score)(args, score_name) for score_name in adni.biomarker_names )
    
def fit_score( args, score_name ):
    print 'Fitting curve to', score_name
    r_file = os.path.join( adni.project_folder, 'src/fitting/vgam_estimate_curves.R' )
    csv_file = os.path.join( adni.project_folder, 'data', score_name.replace( ' ', '_' ) + '.csv' )
    output_file = csv_file.replace( '.csv', '_curves.csv' )
    densities_file = csv_file.replace( '.csv', '_densities.csv' )
    image_file = csv_file.replace( '.csv', '_curves.pdf' )
    stdout_file = csv_file.replace( '.csv', '_stdout.Rout' )

    command = "R CMD BATCH \"--args input_file='%s' output_file='%s' degrees_of_freedom=%i save_plot=1 plot_file='%s' densities_file='%s'\" %s %s" \
              % (csv_file, output_file, int(args.degrees_of_freedom), image_file, densities_file, r_file, stdout_file)
        
    subprocess.call(command, shell=True)
    
def main():
    # parse command line options
    parser = argparse.ArgumentParser(description='Generate reference curves and reports for BAMBI.')
    parser.add_argument( '-n', '--nr_threads', dest = 'nr_threads', type=int, default=4, help='number of threads' )
    parser.add_argument( '-d', '--degrees-of-freedom', dest='degrees_of_freedom', type=int, default=2, help='degrees of freedom for the LMS method')
    parser.add_argument( '--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    #generate_csv_files()
    fit_data( args )

if __name__ == "__main__":
    main()
