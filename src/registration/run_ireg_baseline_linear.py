#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import joblib as jl
import common.adni_tools as adni
import ireg_linear

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( '-n', '--nr_threads', dest = 'nr_threads', type=int, default = 1 )
a = parser.parse_args()

rreg_params = os.path.join( adni.param_folder, 'params-ireg-rigid.txt' )
areg_params = os.path.join( adni.param_folder, 'params-ireg-affine.txt' )

data_folder = os.path.join( adni.data_folder, a.study )
mask_folder = os.path.join( data_folder, 'native/masks_brain' )

output_folder = adni.make_dir( data_folder, 'baseline_linear' )
output_folder_img = adni.make_dir( output_folder, 'images_unstripped' )
output_folder_dof = adni.make_dir( output_folder, 'dof' )

baseline_folder = os.path.join( data_folder, 'native/images_unstripped' )
followup_folder = os.path.join( data_folder, 'native/images_unstripped' )
baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, a.study, a.viscode )

def run( index ):
    target = baseline_files[index]
    source = followup_files[index]
    source_base = os.path.basename( source )
    
    out_dof = os.path.join( output_folder_dof, source_base.replace('.nii.gz', '.dof.gz') )
    out_warped = os.path.join( output_folder_img, source_base )

    ireg_linear.run( source, target, None, out_dof, rreg_params, areg_params, out_warped )

print 'Found', len( baseline_files ), 'image pairs...'
jl.Parallel( n_jobs=a.nr_threads )( jl.delayed(run)(i) for i in range(len(baseline_files)) )
   
