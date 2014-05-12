#! /usr/bin/env python
# print __doc__

import argparse
import os.path
from subprocess import call
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'field_strength', type=str,  help='the field strength, usually 1.5 for ADNI1 and 3 otherwise' )
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( '-r', '--rid', type=str, default=None )
parser.add_argument( '-s', '--spacing', dest='sx', type=str, default='10' )
a = parser.parse_args()
    
rview = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/rview'

data_folder = os.path.join( adni.data_folder, a.study )
baseline_folder = os.path.join( data_folder, 'native/images' )
if a.trans == 'lin':
    followup_folder = os.path.join( data_folder, 'native/images' )
    dof_folder = os.path.join( data_folder, 'baseline_linear/dof' )
else:
    followup_folder = os.path.join( data_folder, 'baseline_linear/images' )
    dof_folder = os.path.join( data_folder, 'baseline_' + a.trans + '_' + a.sx + 'mm_after_linear/dof' )

baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, a.study, a.field_strength, a.viscode )

print 'Found ' + str(len( baseline_files )) + ' images:'
for i in range( len( baseline_files ) ):
    target = baseline_files[i]
    source = followup_files[i]
    source_base = os.path.basename( source )
    dof = os.path.join( dof_folder, source_base ).replace('.nii.gz', '.dof.gz')

    if a.rid == None or source.find( '_S_' + a.rid ) > 0:
        print '--------------------'
        print 'Target: ' + target
        print 'Source: ' + source
        print 'DOF:    ' + dof
        
        call([ rview, target, source, dof, '-res', '1.5', '-mix' ])
