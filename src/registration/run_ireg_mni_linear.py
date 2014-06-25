#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import joblib as jl
import common.adni_tools as adni
import ireg_linear

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( '-n', '--nr_threads', dest = 'nr_threads', type=int, default = 1 )
a = parser.parse_args()

rreg_params = os.path.join( adni.param_folder, 'params-ireg-rigid.txt' )
areg_params = os.path.join( adni.param_folder, 'params-ireg-affine.txt' )

data_folder = os.path.join( adni.data_folder, a.study )
baseline_folder = os.path.join( data_folder, 'native/images' )
baseline_files = adni.get_baseline( baseline_folder, a.study )

output_folder = adni.make_dir( data_folder, 'MNI152_linear' )
output_folder_img = adni.make_dir( output_folder, 'images' )
output_folder_dof = adni.make_dir( output_folder, 'dof' )

def run( index ):
    source = baseline_files[index] 
    source_base = os.path.basename( source )
    
    out_dof = os.path.join( output_folder_dof, source_base.replace('.nii.gz', '.dof.gz') )
    out_warped = os.path.join( output_folder_img, source_base )
    
    ireg_linear.run( source, adni.mni_atlas, None, out_dof, rreg_params, areg_params, out_warped )

print 'Found', len( baseline_files ), 'images...'
jl.Parallel( n_jobs=a.nr_threads )( jl.delayed(run)(i) for i in range(len(baseline_files)) )
