#!/usr/bin/env python

import argparse
import os.path
import csv
import re
import numpy as np
import common.adni_tools as adni
import subprocess

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
                measurements[rid][viscode].update( {'time' : scantime - time_convert} )

    #
    # Remove non-valid rids        
    measurements = { rid : value for rid, value in measurements.items() if rid in valid_rids }
    
    return measurements

def generate_csv():
    measurements = collect_measurements()
    for score_name in adni.biomarker_names:
        print 'Generating output CSV for', score_name
        
        csv_file = os.path.join( adni.project_folder, 'data', score_name.replace( ' ', '_' ) + '.csv' )
        writer = csv.writer( open( csv_file, 'wb' ), delimiter=',' )
        writer.writerow(['rid', 'scan_age', 'volume'])

        for rid, rid_data in measurements.items():
            for _, scan_data in rid_data.items():
                writer.writerow( [rid, scan_data['time'], scan_data[score_name]] )


def generate_reference_data( args ):
    commands = []
    
    for score_name in adni.biomarker_names:
        print 'Fitting curve to', score_name
        
        csv_file = os.path.join( adni.project_folder, 'data', score_name.replace( ' ', '_' ) + '.csv' )
        output_file = csv_file.replace( '.csv', '_curves.csv' )
        image_file = csv_file.replace( '.csv', '_curves.pdf' )
        stdout_file = csv_file.replace( '.csv', '_stdout.Rout' )

        commands += ["R CMD BATCH \"--args input_file='%s' output_file='%s' degrees_of_freedom=%i save_plot=1 plot_filename='%s'\" ref_curve_generator.R %s" 
                     % (csv_file, output_file, int(args.degrees_of_freedom), image_file, stdout_file)]

    for command in commands:
        print command
        return subprocess.call(command, shell=True)
    

def main():
    # parse command line options
    parser = argparse.ArgumentParser(description='Generate reference curves and reports for BAMBI.')
    parser.add_argument('--measurements-path', action='store', dest='measurements_path', help='Path to where the measurements files are collected.')
    parser.add_argument('--output-file', action='store', dest='output_file', help='Pattern of the outputfile.')
    parser.add_argument('--output-curve-dir', action='store', dest='output_curve_dir', help='Directory where the R curves are saved')
    parser.add_argument('--info-file', action='store', dest='info_file', help='Path to the file containing the age and gender information of the subjects.')
    parser.add_argument('--generate-total-report', action='store_true', dest='generate_total_report', default=False, help='Create a single file with all information.')
    parser.add_argument('--degrees-of-freedom', action='store', dest='degrees_of_freedom', type=int, default=2, help='Degrees of freedom for the LMS method.')
    parser.add_argument('--no-separate-csv', action='store_true', dest='no_separate_csv', default=False, help='Do not generate separate csv files.')
    parser.add_argument('--no-percentile-curves', action='store_true', dest='no_percentile_curves', default=False, help='Do not generate percentile curves.')
    parser.add_argument('--generate-single-curve-id', action='store', dest='generate_single_curve_id', help='generate percentile curves for a single output file. Input in the form (abs|picv|ptbv)_(female|male)_Region(x)_(tissue) i.e. abs_male_Region5_wm')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    #generate_csv()
    generate_reference_data( args )

if __name__ == "__main__":
    main()
