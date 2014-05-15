#! /usr/bin/env python
# print __doc__

import argparse
import os.path
from subprocess import call
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
a = parser.parse_args()

execMask = 'padding'

data_folder = os.path.join( adni.data_folder, a.study )
mask_folder = os.path.join( data_folder, 'native/masks_brain' )

baseline_folder = os.path.join( data_folder, 'native/images_unstripped' )
followup_folder = os.path.join( data_folder, 'baseline_linear/images_unstripped' )
baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, a.study, a.viscode )

output_folder = adni.make_dir( data_folder, 'baseline_linear/images' )

print 'Found ' + str(len( baseline_files )) + ' images...'
for i in range( len( baseline_files ) ):
    baseline = baseline_files[i]
    baseline_base = os.path.basename( baseline )
    followup = followup_files[i]
    followup_base = os.path.basename( followup )
    
    mask = os.path.join( mask_folder, baseline_base )
    mask = adni.check_mask( mask )
    if mask != None:
        mask_base = os.path.basename( mask )
        
        out = os.path.join( output_folder, followup_base )
        
        if os.path.isfile( out ):
            print 'Image already exists: ' + out
        elif mask == None or not os.path.isfile( mask ):
            print 'No mask found for: ' + out
        else:
            print '--------------------'
            print 'Image:  ' + followup
            print 'Mask:   ' + mask
            print 'Output: ' + out
            call([ execMask, followup, mask, out, '1', '0', '-invert'])
