#! /usr/bin/env python
# print __doc__

import sys
import os.path
import threading
import common.adni_tools as adni
import ireg_nonlinear

if len( sys.argv ) < 6:
    print 'Usage: run_ireg_mni_nonlinear_after_linear.py <threads> <study> <field_strength> <transformation> <spacing>'
    exit()
  
nr_threads = int( sys.argv[1] )
study = sys.argv[2]
fs = sys.argv[3]
trans = sys.argv[4]
sx = sys.argv[5]
viscode = 'bl'

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', study )

target_mni = '/vol/medic01/users/aschmidt/projects/Data/MNI152-Template/MNI152_T1_1mm_brain.nii'
ireg_params = '/vol/biomedic/users/aschmidt/ADNI/scripts/registration/params-ireg-' + trans + '-' + sx + 'mm.txt'

if viscode == 'bl':
    baseline_folder = os.path.join( data_folder, 'MNI152_linear/images' )
    baseline_files = adni.get_baseline( baseline_folder, study, fs, 'bl' )
    output_folder = adni.make_dir( data_folder, 'MNI152_' + trans + '_' + sx + 'mm_after_linear' )
elif viscode == 'm12' or viscode == 'm24':
    baseline_folder = os.path.join( data_folder, 'MNI152_linear_via_baseline/images' )
    baseline_files = adni.get_baseline( baseline_folder, study, fs, 'm24' )
    output_folder = adni.make_dir( data_folder, 'MNI152_' + trans + '_' + sx + 'mm_after_linear_via_baseline' )

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
     
    if thread_ctr == nr_threads:
        for t in threads:
            t.join()
        threads = []
        thread_ctr = 0
