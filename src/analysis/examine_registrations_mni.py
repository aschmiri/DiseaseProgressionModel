#! /usr/bin/env python
# print __doc__

import sys
import os.path
from subprocess import call
import common.adni_tools as adni

if len( sys.argv ) < 5:
    print 'Usage: examine_registrations_mni.py <study> <field_strength> <transformation> <spacing> [<subject>]'
    exit()
  
study = sys.argv[1]
fs = sys.argv[2]
trans = sys.argv[3]
sx = sys.argv[4]

subject = None
if len( sys.argv ) == 6:
    subject = sys.argv[5]

with_dof = False

rview = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/rview'
target_mni = '/vol/medic01/users/aschmidt/projects/Data/MNI152-Template/MNI152_T1_1mm_brain_float.nii'

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', study )

if with_dof:
    baseline_folder = os.path.join( data_folder, 'MNI152_linear/images' )
else:
    baseline_folder = os.path.join( data_folder, 'MNI152_' + trans + '_' + sx + 'mm_after_linear/images' )
dof_folder = os.path.join( data_folder, 'MNI152_' + trans + '_' + sx + 'mm_after_linear/dof' )

baseline_files = adni.get_baseline( baseline_folder, study, fs )
baseline_files, dof_files = adni.find_images_with_dof( baseline_files, dof_folder )

print 'Found ' + str(len( baseline_files )) + ' images:'
for i in range( len( baseline_files ) ):
    source = baseline_files[i]
    dof = dof_files[i]
    if subject == None or source.find( subject ) > 0:
        print '--------------------'
        print 'Source: ' + source
        print 'DOF:    ' + dof
        
        if with_dof:
            call([ rview, target_mni, source, dof, '-res', '1.5', '-mix' ])
        else:
            call([ rview, target_mni, source, '-res', '1.5', '-mix' ])
        
        
# ADNI_127_S_1032 3
# ADNI_018_S_0103 3
# ADNI_023_S_1046 2
#
#
#
#
#
#
# 
