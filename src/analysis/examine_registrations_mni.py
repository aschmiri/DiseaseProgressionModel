#! /usr/bin/env python
# print __doc__

import argparse
import os.path
from subprocess import call
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'field_strength', type=str,  help='the field strength, usually 1.5 for ADNI1 and 3 otherwise' )
parser.add_argument( 'trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( '-d', '--dof', type=bool, target='with_dof', action='store_true', default=False, help='show the dof' )
parser.add_argument( '-r', '--rid', type=int, default=None )
parser.add_argument( '-s', '--spacing', dest='sx', type=str, default='10' )
a = parser.parse_args()

rview = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/rview'
target_mni = os.path.join( adni.mni_folder, 'MNI152_T1_1mm_brain.nii' )

data_folder = os.path.join( adni.data_folder, a.study )
if a.with_dof:
    baseline_folder = os.path.join( data_folder, 'MNI152_linear/images' )
else:
    baseline_folder = os.path.join( data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear/images' )
dof_folder = os.path.join( data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear/dof' )

baseline_files = adni.get_baseline( baseline_folder, a.study, a.fs )
baseline_files, dof_files = adni.find_images_with_dof( baseline_files, dof_folder )

print 'Found ' + str(len( baseline_files )) + ' images:'
for i in range( len( baseline_files ) ):
    source = baseline_files[i]
    dof = dof_files[i]
    if a.rid == None or source.find( '_S_' + str( a.rid ) ) > 0:
        print '--------------------'
        print 'Source: ' + source
        print 'DOF:    ' + dof
        
        if a.with_dof:
            call([ rview, target_mni, source, dof, '-res', '1.5', '-mix' ])
        else:
            call([ rview, target_mni, source, '-res', '1.5', '-mix' ])
