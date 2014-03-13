#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import threading
import common.adni_tools as adni
import ireg_nonlinear

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'field_strength', type=str,  help='the field strength, usually 1.5 for ADNI1 and 3 otherwise' )
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( '-n', '--nr_threads', dest = 'nr_threads', type=int, default = 1 )
parser.add_argument( '-s', '--spacing', dest = 'sx', type=str, default = '10' )
a = parser.parse_args()

ireg_params = os.path.join( adni.param_folder, 'params-ireg-' + a.trans + '-' + a.sx + 'mm.txt' )
target_mni = os.path.join( adni.mni_folder, 'MNI152_T1_1mm_brain.nii' )

data_folder = os.path.join( adni.data_folder, a.study )
if a.viscode == 'bl':
    baseline_folder = os.path.join( data_folder, 'MNI152_linear/images' )
    baseline_files = adni.get_baseline( baseline_folder, a.study, a.field_strength, 'bl' )
    output_folder = adni.make_dir( data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear' )
elif a.viscode in ['m06', 'm12', 'm24']:
    baseline_folder = os.path.join( data_folder, 'MNI152_linear_via_baseline/images' )
    baseline_files = adni.get_baseline( baseline_folder, a.study, a.field_strength, 'm24' )
    output_folder = adni.make_dir( data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear_via_baseline' )

output_folder_img = adni.make_dir( output_folder, 'images' )
output_folder_dof = adni.make_dir( output_folder, 'dof' )

class RegistrationThread(threading.Thread):
    def __init__(self, index):
        threading.Thread.__init__(self)
        self.index = index
    def run(self):
        source = baseline_files[self.index] 
        source_base = os.path.basename( source )
        
        out_dof = os.path.join( output_folder_dof, source_base.replace('.nii.gz', '.dof.gz') )
        out_warped = os.path.join( output_folder_img, source_base )
        
        ireg_nonlinear.run( source, target_mni, 'none', out_dof, ireg_params, out_warped )

print 'Found ' + str(len( baseline_files )) + ' images...'
thread_ctr = 0
threads = []
for i in range( len( baseline_files ) ):
    thread = RegistrationThread(i)
    thread.start()
    threads.append(thread)
    thread_ctr += 1
     
    if thread_ctr == a.nr_threads:
        for t in threads:
            t.join()
        threads = []
        thread_ctr = 0
