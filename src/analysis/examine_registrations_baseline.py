#! /usr/bin/env python
# print __doc__

import sys
import os.path
from subprocess import call
import common.adni_tools as adni

if len( sys.argv ) < 5:
    print 'Usage: examine_registrations_baseline.py <study> <field_strength> <transformation> <spacing> [<subject>]'
    exit()
  
study = sys.argv[1]
fs = sys.argv[2]
trans = sys.argv[3]
sx = sys.argv[4]

subject = None
if len( sys.argv ) == 6:
    subject = sys.argv[5]
    
rview = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/rview'

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', study )

#baseline_folder = os.path.join( data_folder, 'native/images_unstripped' )
#followup_folder = os.path.join( data_folder, 'baseline_linear/images_unstripped' )
baseline_folder = os.path.join( data_folder, 'native/images' )
followup_folder = os.path.join( data_folder, 'baseline_linear/images' )
dof_folder = os.path.join( data_folder, 'baseline_' + trans + '_' + sx + 'mm_after_linear/dof' )

baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, study, fs )

print 'Found ' + str(len( baseline_files )) + ' images:'
for i in range( len( baseline_files ) ):
    target = baseline_files[i]
    source = followup_files[i]
    source_base = os.path.basename( source )
    dof = os.path.join( dof_folder, source_base )

    if subject == None or source.find( subject ) > 0:
        print '--------------------'
        print 'Target: ' + target
        print 'Source: ' + source
        print 'DOF:    ' + dof
        
        call([ rview, target, source, dof, '-res', '1.5', '-mix' ])
