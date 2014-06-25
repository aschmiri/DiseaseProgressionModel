#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import joblib as jl
import common.adni_tools as adni
import common.atlas_tools as at
import ireg_nonlinear

parser = argparse.ArgumentParser()
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...' )
parser.add_argument( 'trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( 'state', type=float, help='the state for which relevant images should be registered' )
parser.add_argument( '-n', '--nr_threads', type=int, default=1 )
parser.add_argument( '-i', '--iteration', type=int, default=1 )
parser.add_argument( '-r', '--required_subjects', type=int, default=20 )
parser.add_argument( '-s', '--spacing', type=str, default='10' )
parser.add_argument( '-a', '--age_regression', action='store_true', default=False, help='use age regression' )
parser.add_argument( '--save_image', action='store_true', help='save the warped image' )
a = parser.parse_args()

ireg_params = os.path.join( adni.param_folder, 'params-ireg-' + a.trans + '-' + a.spacing + 'mm.txt' )

atlas_folder = os.path.join( adni.project_folder, 'atlas/model_' + str(a.iteration) )
datafile = os.path.join( atlas_folder, 'data_' + a.trans +'_' + a.viscode + '_' + a.diagnosis + '.csv' )

rids, _, _, states, images = at.read_datafile( datafile, a.diagnosis, age_regression = a.age_regression )

sigma, weights, indices = at.adaptive_kernel_regression( states, a.state, required_subjects = a.required_subjects )

selected_rids = rids[indices]
selected_images = images[indices]
selected_weights = weights[indices]

output_folder = adni.make_dir( adni.data_folder, 'ADNI/MNI152_intra_' + a.trans + '_' + a.spacing + 'mm' )
output_folder_dof = adni.make_dir( output_folder, 'dof' )
if a.save_image:
    output_folder_img = adni.make_dir( output_folder, 'images' )
    
def run( index_targ, index_srce ):
    target = selected_images[index_targ]
    source = selected_images[index_srce]
    target_rid = selected_rids[index_targ]
    source_rid = selected_rids[index_srce]
    out_basename = str(source_rid) + '_to_' + str(target_rid)

    out_dof = os.path.join( output_folder_dof, out_basename + '.dof.gz' )
    if a.save_image:
        out_warped = os.path.join( output_folder_img, out_basename + '.nii.gz' )
    else:
        out_warped = None
            
    ireg_nonlinear.run( source, target, None, out_dof, ireg_params, out_warped )

print 'Found', len( selected_images ), 'relevant images for state', a.state, '...'
for j in range(len(selected_images)):
    jl.Parallel( n_jobs=a.nr_threads )( 
        jl.delayed(run)(i,j) for i in range(len(selected_images)) )
