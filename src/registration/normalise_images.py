#! /usr/bin/env python
# print __doc__

import argparse
import os.path
from subprocess import call
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'study' )
a = parser.parse_args()

#reg_algorithm = 'MNI152_svffd_10mm_after_linear'
reg_algorithm = 'MNI152_linear'

execNormalise = 'normalize_percentiles_apply'
percentiles = '/vol/biomedic/users/reg09/posdoc/manifold_alignment/percentiles_10.csv'

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
reg_folder = os.path.join( base_folder, 'data', a.study, reg_algorithm )
image_folder = os.path.join( reg_folder, 'images' )
if not os.path.isdir( image_folder ):
    print 'Input folder does not exist.'
    exit()

norm_folder = adni.make_dir( reg_folder, 'images_normalised' )

for infile in os.listdir( image_folder ):
    if infile.endswith( '.nii.gz' ):
        image_in = os.path.join( image_folder, infile )
        image_out = os.path.join( norm_folder, infile )
        
        if os.path.exists( image_out ):
            print 'Image already exists: ' + image_out
        else:
            print '--------------------'
            print 'Starting: normalize_percentiles_apply'
            print 'Input:  ' + image_in
            print 'Output: ' + image_out
            call([ execNormalise, percentiles, image_in, norm_folder + '/' ])
        