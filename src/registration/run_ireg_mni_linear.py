#! /usr/bin/env python
# print __doc__

import sys
import os.path
import threading
import common.adni_tools as adni
import ireg_linear

if len( sys.argv ) < 4:
    print 'Usage: run_ireg_mni_linear.py <threads> <study> <field_strength>'
    exit()
  
nr_threads = int( sys.argv[1] )
study = sys.argv[2]
fs = sys.argv[3]

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', study )

target_mni = '/vol/medic01/users/aschmidt/projects/Data/MNI152-Template/MNI152_T1_1mm_brain.nii'
rreg_params = '/vol/biomedic/users/aschmidt/ADNI/scripts/registration/params-ireg-rigid.txt'
areg_params = '/vol/biomedic/users/aschmidt/ADNI/scripts/registration/params-ireg-affine.txt'

baseline_folder = os.path.join( data_folder, 'native/images' )
baseline_files = adni.get_baseline( baseline_folder, study, fs )

output_folder = adni.make_dir( data_folder, 'MNI152_linear' )
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
        
        ireg_linear.run( source, target_mni, 'none', out_dof, rreg_params, areg_params, out_warped )

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
