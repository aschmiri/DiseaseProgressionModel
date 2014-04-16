#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import threading
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
    
class RegistrationThread(threading.Thread):
    def __init__(self, index_targ, index_srce):
        threading.Thread.__init__(self)
        self.index_targ = index_targ
        self.index_srce = index_srce
    def run(self):
        target = selected_images[self.index_targ]
        source = selected_images[self.index_srce]
        target_rid = selected_rids[self.index_targ]
        source_rid = selected_rids[self.index_srce]
        out_basename = str(source_rid) + '_to_' + str(target_rid)

        out_dof = os.path.join( output_folder_dof, out_basename + '.dof.gz' )
        if a.save_image:
            out_warped = os.path.join( output_folder_img, out_basename + '.nii.gz' )
        else:
            out_warped = None
                
        ireg_nonlinear.run( source, target, None, out_dof, ireg_params, out_warped )

print 'Found ' + str(len( selected_images )) + ' relevant images for state ' + str(a.state) + '...'
thread_ctr = 0
threads = []
for i in range( len( selected_images ) ):
    for j in range( len( selected_images ) ):
        if i != j:
            thread = RegistrationThread( i, j )
            thread.start()
            threads.append(thread)
            thread_ctr += 1
             
            if thread_ctr == a.nr_threads:
                for t in threads:
                    t.join()
                threads = []
                thread_ctr = 0
