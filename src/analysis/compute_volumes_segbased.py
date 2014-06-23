#! /usr/bin/env python
# print __doc__

import os.path
import argparse
import csv
import sqlite3
from subprocess import check_output
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( '-m', '--min_scans', dest = 'min_scans', type=int, default=6, help='the minimal number of scans per subject' )
a = parser.parse_args()

exec_volumes = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_compute_volume'
exec_factors = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/asr_affineScaling'

data_file = os.path.join( adni.project_folder, 'lists/query_ADNI.csv' )
out_file = os.path.join( adni.project_folder, 'lists/volumes_segbased.csv' )
folder = adni.data_folder + '/ALL/native/images_unstripped'

folder_lin_mni = adni.data_folder + '/ALL/MNI152_linear/dof'
folder_lin_bl  = adni.data_folder + '/ALL/baseline_linear/dof'

#
# Read csv file with metadata
#
con = sqlite3.connect( os.path.join( adni.project_folder, 'lists', 'adni.db' ) )
con.row_factory = sqlite3.Row

cur = con.cursor()
cur.execute( "SELECT rid FROM adnimerge WHERE viscode = 'bl'" )
rows = cur.fetchall()
rids = [row['rid'] for row in rows]

with open( out_file, 'wb' ) as csvfile:
    writer = csv.writer( csvfile, delimiter=',' )
    writer.writerow( ['RID', 'VISCODE', 'DX.scan', 'FactorMNI','FactorBL'] + adni.volume_names )
    
    for rid in rids:
        cur.execute( "SELECT iid, viscode, study, diagnosis, filename FROM adnimerge WHERE rid = " + str(rid) )
        scans = cur.fetchall()
        if len(scans) >= a.min_scans:
            for scan in scans:
                viscode = scan['viscode']
                iid = scan['iid']
                study = scan['study']
                filename = scan['filename']
                diagnosis = scan['diagnosis']

                if viscode == 'bl':
                    seg = adni.find_file( os.path.join( adni.data_folder, study, 'native/seg_138regions_baseline', 'EM-' + filename ) )
                    mni_dof = os.path.join( folder_lin_mni.replace( 'ALL', study ), filename.replace( '.nii.gz', '.dof.gz' ) )
                    if not os.path.isfile( mni_dof ): 
                        print 'ERROR: File not found:', mni_dof
                        mni_factor = -1
                    else:
                        mni_factor = float( check_output([exec_factors, mni_dof]) )
                    bl_factor = 1
                        
                else:
                    # Get baseline file
                    cur.execute( "SELECT study, filename FROM adnimerge WHERE rid = " + str(rid) + " AND viscode = 'bl'" )
                    bl_scans = cur.fetchall()
                    if len(bl_scans) != 1:
                        print 'WARNING: wrong number of baseline files found for', rid
                    else:
                        bl_study = bl_scans[0]['study']
                        bl_base = bl_scans[0]['filename']

                    seg = adni.find_file( os.path.join( adni.data_folder, study, 'native/seg_138regions_followup_warped', filename ) )
                    mni_dof = os.path.join( folder_lin_mni.replace( 'ALL', bl_study ), bl_base.replace( '.nii.gz', '.dof.gz' ) )
                    bl_dof = os.path.join( folder_lin_bl.replace( 'ALL', study ), filename.replace( '.nii.gz', '.dof.gz' ) )
                    if not os.path.isfile( mni_dof ): 
                        print 'ERROR: File not found:', mni_dof
                        mni_factor = -1
                    else:
                        mni_factor = float( check_output([exec_factors, mni_dof]) )
                    
                    if not os.path.isfile( bl_dof ):
                        print 'ERROR: File not found:', bl_dof
                        bl_factor = -1
                    else:
                        bl_factor = float( check_output([exec_factors, bl_dof]) )
                        
                # Get volumes of 138 objects
                if seg == None or not os.path.isfile( seg ):
                    print 'ERROR: File not found:', seg
                elif bl_factor == -1: # or mni_factor == -1:
                    print 'ERROR: Factors were not computed correctly'
                else:
                    volumes = check_output([exec_volumes, seg])
                    volumes = [float(vol) for vol in volumes.split(',')]
                    volumes.pop(0)
                    if len(volumes) != 138:
                        print 'ERROR:', len(volumes), 'volumes read for', os.path.basename( seg )
                    else:
                        writer.writerow( [str(rid), viscode, diagnosis, mni_factor, bl_factor] + volumes )
                        csvfile.flush()
                        #print rid, viscode, diagnosis, mni_factor, bl_factor, ':', volumes
