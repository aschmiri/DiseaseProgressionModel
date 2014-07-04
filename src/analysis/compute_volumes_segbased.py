#! /usr/bin/env python
# print __doc__

import os.path
import argparse
import csv
import sqlite3
from subprocess import check_output
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( '-m', '--min_scans', dest='min_scans', type=int, default=6, help='the minimal number of scans per subject' )
parser.add_argument( '-t', '--trans', dest='trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( '-s', '--spacing', dest='sx', type=str, default='5' )
a = parser.parse_args()

exec_volumes = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_compute_volume'
exec_factors = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/asr_affineScaling'

data_file = os.path.join( adni.project_folder, 'lists/query_ADNI.csv' )
out_file = os.path.join( adni.project_folder, 'lists/volumes_segbased_' + a.trans + '_' + a.sx + 'mm__.csv' )
folder = adni.data_folder + '/ALL/native/images_unstripped'

folder_lin_mni = adni.data_folder + '/ALL/MNI152_linear/dof'
folder_lin_bl  = adni.data_folder + '/ALL/baseline_linear/dof'

# Setup DB
con = sqlite3.connect( os.path.join( adni.project_folder, 'lists', 'adni.db' ) )
con.row_factory = sqlite3.Row
cur = con.cursor()

# Get all RIDs
cur.execute( "SELECT rid FROM adnimerge WHERE viscode = 'bl'" )
rows = cur.fetchall()
rids = [row['rid'] for row in rows]

# Estimate volumes
with open( out_file, 'wb' ) as csvfile:
    writer = csv.writer( csvfile, delimiter=',' )
    writer.writerow( ['RID', 'VISCODE', 'DX.scan', 'AGE.scan','ScanDate',
                      'MMSE','CDRSB','ADAS11','ADAS13','FAQ','MOCA',
                      'FactorMNI','FactorBL'] + adni.volume_names )
    
    for rid in rids:
        cur.execute( "SELECT iid, viscode, study, diagnosis, filename, age, scandate \
                      FROM adnimerge WHERE rid = " + str(rid) )
        scans = cur.fetchall()
        if len(scans) >= a.min_scans:
            for scan in scans:
                # Get basic parameters
                viscode   = scan['viscode']
                iid       = scan['iid']
                study     = scan['study']
                filename  = scan['filename']
                diagnosis = scan['diagnosis']
                age       = scan['age']
                scandate  = scan['scandate']
                
                # Get cognitive scores
                cur.execute( "SELECT mmse, cdrsb, adas11, adas13, faq, moca \
                              FROM adnimerge JOIN adnimerge_cog USING (rid, viscode) \
                              WHERE rid = " + str(rid) + " AND viscode = '" + viscode + "'" )
                cog_data = cur.fetchall()
                
                if len(cog_data) != 1:
                    print 'WARNING: cog_data has wrong size:', len(cog_data)
                mmse   = adni.safe_cast( cog_data[0]['mmse'] )
                cdrsb  = adni.safe_cast( cog_data[0]['cdrsb'] )
                adas11 = adni.safe_cast( cog_data[0]['adas11'] )
                adas13 = adni.safe_cast( cog_data[0]['adas13'] )
                faq    = adni.safe_cast( cog_data[0]['faq'] )
                moca   = adni.safe_cast( cog_data[0]['moca'] )

                # Get data for volume computation
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

                    seg = adni.find_file( os.path.join( adni.data_folder, study, 'native/seg_138regions_followup_' + a.trans + '_' + a.sx + 'mm', filename ) )
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
                if seg != None:
                    if bl_factor == -1: # or mni_factor == -1:
                        print 'ERROR: Factors were not computed correctly'
                    else:
                        volumes = check_output([exec_volumes, seg])
                        volumes = [float(vol) for vol in volumes.split(',')]
                        volumes.pop(0)
                        if len(volumes) != 138:
                            print 'ERROR:', len(volumes), 'volumes read for', os.path.basename( seg )
                        else:
                            writer.writerow( [str(rid), viscode, diagnosis, age, scandate,
                                              mmse,cdrsb,adas11,adas13,faq,moca,
                                              mni_factor, bl_factor] + volumes )
                            csvfile.flush()
                            #print rid, viscode, diagnosis, mni_factor, bl_factor, ':', volumes
