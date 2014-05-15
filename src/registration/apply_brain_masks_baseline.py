#! /usr/bin/env python
# print __doc__

import argparse
import os.path
from subprocess import call
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
a = parser.parse_args()

execMask = 'padding'
data_folder = os.path.join( adni.data_folder, a.study )
output_folder = os.path.join( data_folder, 'native/images' )
mask_folder = os.path.join( data_folder, 'native/masks_brain' )
baseline_folder = os.path.join( data_folder, 'native/images_unstripped' )

baseline_files = adni.get_baseline( baseline_folder, a.study )

print 'Found', str(len( baseline_files )), 'images...'
for i in range( len( baseline_files ) ):
    image = baseline_files[i]
    image_base = os.path.basename( image )

    mask = os.path.join( mask_folder, image_base )
    mask = adni.check_mask( mask )
    if mask != None:
        mask_base = os.path.basename( mask )
    
        out = os.path.join( output_folder, image_base )
        
        if os.path.isfile( out ):
            print 'Image already exists:', out
        elif not os.path.isfile( mask ):
            print 'No mask found for: ' + out
        else:
            print '--------------------'
            print 'Image: ', image
            print 'Mask:  ', mask
            print 'Output:', out
             
            call([ execMask, image, mask, out, '1', '0', '-invert'])
