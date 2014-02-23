#! /usr/bin/env python
# print __doc__

import os.path
import argparse
from subprocess import check_output
import csv
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...' )
parser.add_argument( '-i', '--iteration', dest='iteration', type=int, default=0 )
parser.add_argument( '--state_min', type=float, default=0 )
parser.add_argument( '--state_max', type=float, default=15 )
parser.add_argument( '--state_stepsize', type=float, default = 0.25 )
a = parser.parse_args()
                 
execEstimate = 'stateestimation'

mask_brain = '/vol/medic01/users/aschmidt/projects/Data/MNI152-Template/MNI152_T1_1mm_brain.nii'

atlas_folder = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_' + str(a.iteration)
data_folder = '/vol/biomedic/users/aschmidt/ADNI'

image_folder_adni1_bl = os.path.join( data_folder, 'data/ADNI1/MNI152_linear/images' )
image_folder_adni2_bl = os.path.join( data_folder, 'data/ADNI2/MNI152_linear/images' )
image_folder_adniG_bl = os.path.join( data_folder, 'data/ADNIGO/MNI152_linear/images' )
datafile_bl = os.path.join( atlas_folder, 'data_' + a.viscode + '_' + a.diagnosis + '.csv' )

image_folder_adni1_fu = os.path.join( data_folder, 'data/ADNI1/MNI152_linear_via_baseline/images' )
image_folder_adni2_fu = os.path.join( data_folder, 'data/ADNI2/MNI152_linear_via_baseline/images' )
image_folder_adniG_fu = os.path.join( data_folder, 'data/ADNIGO/MNI152_linear_via_baseline/images' )
datafile_m12 = os.path.join( atlas_folder, 'data_' + a.viscode + '_' + a.diagnosis + '_testm12.csv' )
datafile_m24 = os.path.join( atlas_folder, 'data_' + a.viscode + '_' + a.diagnosis + '_testm24.csv' )

model_prefix = os.path.join( atlas_folder, 'model_' + a.viscode + '_' + a.diagnosis + '_s' )

def collect_data( data_folder_adni1, data_folder_adni2, data_folder_adniG,
                   viscode, datafile ):
    
    files, rids, diags, ages, mmses = adni.get_all_data( 
        data_folder_adni1, data_folder_adni2, data_folder_adniG, viscode )
    
    if os.path.exists( datafile ):
        print ' File ' + datafile + ' already exists.'
    else:
        with open( datafile, 'wb') as csvfile:
            csv_writer = csv.writer( csvfile, delimiter=',' )
            csv_writer.writerow([ 'RID', 'DX.bl', 'AGE', 'MMSE', 'VIRT_NMI', 'FILE' ])
            for i in range(len( files )):
                print ' Analysing ' + files[i]
                
                # Estimate disease state for image            
                if os.path.exists( files[i] ) and a.diagnosis in diags[i] or diags[i] in a.diagnosis:
                    virtual_nmi_brain = check_output([execEstimate, files[i], '-mask', mask_brain,
                                                      '-min', str(a.state_min), '-max', str(a.state_max), '-stepsize', str(a.state_stepsize),
                                                      '-prefix', model_prefix, '-silent' ])
                    csv_writer.writerow([rids[i], diags[i], ages[i], mmses[i], virtual_nmi_brain, files[i]])
                    print virtual_nmi_brain

collect_data( image_folder_adni1_bl, image_folder_adni2_bl, image_folder_adniG_bl, 'bl', datafile_bl )
#collect_data( image_folder_adni1_fu, image_folder_adni2_fu, image_folder_adniG_fu, 'm24', datafile_m24 )
#collect_data( image_folder_adni1_fu, image_folder_adni2_fu, image_folder_adniG_fu, 'm12', datafile_m12 )
