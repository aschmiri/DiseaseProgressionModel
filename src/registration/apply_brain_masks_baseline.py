#! /usr/bin/env python
# print __doc__

import sys
import os.path
from subprocess import call
import common.adni_tools as adni

if len( sys.argv ) < 3:
    print 'Usage: apply_brain_masks_baseline.py <study> <field_strength>'
    exit()
  
study = sys.argv[1]
fs = sys.argv[2]

execMask = 'padding'

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', study )
output_folder = os.path.join( data_folder, 'native/images' )
mask_folder = os.path.join( data_folder, 'native/masks_brain' )
baseline_folder = os.path.join( data_folder, 'native/images_unstripped' )

baseline_files = adni.get_baseline( baseline_folder, study, fs )

print 'Found ' + str(len( baseline_files )) + ' images...'
for i in range( len( baseline_files ) ):
    image = baseline_files[i]
    image_base = os.path.basename( image )

    mask = os.path.join( mask_folder, image_base )
    mask = adni.check_mask( mask )
    mask_base = os.path.basename( mask )

    out = os.path.join( output_folder, image_base )
    
    if os.path.isfile( out ):
        print 'Image already exists: ' + out
    else:
        print '--------------------'
        print 'Image:  ' + image
        print 'Mask:   ' + mask
        print 'Output: ' + out
         
        call([ execMask, image, mask, out, '1', '0', '-invert'])
