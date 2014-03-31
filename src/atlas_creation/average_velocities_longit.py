#! /usr/bin/env python
# print __doc__

import os.path
import argparse
from subprocess import call
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...' )
a = parser.parse_args()

execAverage = 'ffdaverage'

dof_folder_adni1 = os.path.join( adni.data_folder, 'ADNI1', 'MNI152_svffd_10mm_followup_to_baseline/dof' )
dof_folder_adni2 = os.path.join( adni.data_folder, 'ADNI2', 'MNI152_svffd_10mm_followup_to_baseline/dof' )

velocities_m12, rids_m12 = adni.get_baseline_transformations_and_rids( dof_folder_adni1, dof_folder_adni2, 'm12', a.diagnosis )
velocities_m24, rids_m24 = adni.get_baseline_transformations_and_rids( dof_folder_adni1, dof_folder_adni2, 'm24', a.diagnosis )

atlas_folder = os.path.join( adni.project_folder, 'atlas/model_0_longit' )
out_average_velo_m12 = os.path.join( atlas_folder, 'velo_m12_' + a.diagnosis + '.dof.gz' )  
out_average_velo_m24 = os.path.join( atlas_folder, 'velo_m24_' + a.diagnosis + '.dof.gz' )

selected_velocities_m12 = []
selected_velocities_m24 = []

for i in range( len(rids_m12) ):
    rid_m12 = rids_m12[i]
    for j in range( len(rids_m24) ):
        rid_m24 = rids_m24[j]
        if rid_m24 == rid_m12:
            selected_velocities_m12.append( velocities_m12[i] )
            selected_velocities_m24.append( velocities_m24[j] )

def average_velovities( velocities, out_average_velo ):
    print 'Found ' + str(len( velocities )) + ' velocities in total...'
    if os.path.isfile( out_average_velo ):
        print 'Image already exists: ' + out_average_velo
    else:
        print '--------------------'
        print 'Output velocities: ' + out_average_velo
          
        call([ execAverage, out_average_velo] + velocities )

average_velovities( selected_velocities_m12, out_average_velo_m12 )
average_velovities( selected_velocities_m24, out_average_velo_m24 )
