#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import joblib as jl
import common.adni_tools as adni
import ireg_nonlinear

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( '-n', '--nr_threads', dest = 'nr_threads', type=int, default = 1 )
parser.add_argument( '-s', '--spacing', dest = 'sx', type=str, default = '10' )
a = parser.parse_args()

ireg_params = os.path.join( adni.param_folder, 'params-ireg-' + a.trans + '-' + a.sx + 'mm.txt' )

data_folder = os.path.join( adni.data_folder, a.study )
baseline_folder = os.path.join( data_folder, 'native/images_unstripped' )
followup_folder = os.path.join( data_folder, 'baseline_linear/images_unstripped' )
baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, a.study, a.viscode )

output_folder = adni.make_dir( data_folder, 'baseline_' + a.trans + '_' + a.sx + 'mm_after_linear' )
output_folder_img = adni.make_dir( output_folder, 'images_unstripped' )
output_folder_dof = adni.make_dir( output_folder, 'dof' )

def run( index ):
    target = baseline_files[index]
    source = followup_files[index]
    source_base = os.path.basename( source )
    
    out_dof = os.path.join( output_folder_dof, source_base.replace('.nii.gz', '.dof.gz') )
    out_warped = os.path.join( output_folder_img, source_base )
    
    ireg_nonlinear.run( source, target, None, out_dof, ireg_params, out_warped )

print 'Found', len( baseline_files ), 'image pairs...'
jl.Parallel( n_jobs=a.nr_threads )( jl.delayed(run)(i) for i in range(len(baseline_files)) )
