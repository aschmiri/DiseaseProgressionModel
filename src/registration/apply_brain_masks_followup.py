#! /usr/bin/env python
# print __doc__

import sys
import os.path
from subprocess import call
import common.adni_tools as adni

if len( sys.argv ) < 4:
    print 'Usage: apply_brain_masks_followup.py <study> <field_strength> <viscode>'
    exit()
  
study = sys.argv[1]
fs = sys.argv[2]
viscode = sys.argv[3]

execMask = 'padding'

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', study )

mask_folder = os.path.join( data_folder, 'native/masks_brain' )

baseline_folder = os.path.join( data_folder, 'native/images_unstripped' )
followup_folder = os.path.join( data_folder, 'baseline_linear/images_unstripped' )
baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, study, fs, viscode )

output_folder = adni.make_dir( data_folder, 'baseline_linear/images' )

print 'Found ' + str(len( baseline_files )) + ' images...'
for i in range( len( baseline_files ) ):
    baseline = baseline_files[i]
    baseline_base = os.path.basename( baseline )
    mask = os.path.join( mask_folder, baseline_base )
    mask = adni.check_mask( mask )
    
    mask_base = os.path.basename( mask )
    followup = followup_files[i]
    followup_base = os.path.basename( followup )
    out = os.path.join( output_folder, followup_base )
    
    if os.path.isfile( out ):
        print 'Image already exists: ' + out
    elif mask == 'none':
        print 'No mask found for: ' + out
    else:
        print '--------------------'
        print 'Image:  ' + followup
        print 'Mask:   ' + mask
        print 'Output: ' + out
        call([ execMask, followup, mask, out, '1', '0', '-invert'])
