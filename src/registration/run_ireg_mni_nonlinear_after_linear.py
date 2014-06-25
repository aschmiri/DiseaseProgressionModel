#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import joblib as jl
import common.adni_tools as adni
import ireg_nonlinear
import re

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( '-n', '--nr_threads', dest = 'nr_threads', type=int, default = 1 )
parser.add_argument( '-s', '--spacing', dest = 'sx', type=str, default = '10' )
a = parser.parse_args()

ireg_params = os.path.join( adni.param_folder, 'params-ireg-' + a.trans + '-' + a.sx + 'mm.txt' )

data_folder = os.path.join( adni.data_folder, a.study )
if a.viscode == 'bl':
    image_folder = os.path.join( data_folder, 'MNI152_linear/images' )
    image_files = adni.get_baseline( image_folder, a.study )
    output_folder = adni.make_dir( data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear' )
elif a.viscode == 'fu' or re.match('m[0-9][0-9]', a.viscode):
    image_folder = os.path.join( data_folder, 'MNI152_linear_via_baseline/images' )
    image_files = adni.get_followup( image_folder, a.study, a.viscode )
    output_folder = adni.make_dir( data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear_via_baseline' )

output_folder_img = adni.make_dir( output_folder, 'images' )
output_folder_dof = adni.make_dir( output_folder, 'dof' )

def run( index ):
    source = image_files[index] 
    source_base = os.path.basename( source )
    
    out_dof = os.path.join( output_folder_dof, source_base.replace('.nii.gz', '.dof.gz') )
    out_warped = os.path.join( output_folder_img, source_base )
    
    ireg_nonlinear.run( source, adni.mni_atlas, None, out_dof, ireg_params, out_warped )

print 'Found', len( image_files ), 'images...'
jl.Parallel( n_jobs=a.nr_threads )( jl.delayed(run)(i) for i in range(len(image_files)) )
