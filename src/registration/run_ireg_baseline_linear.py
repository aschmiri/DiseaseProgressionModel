#! /usr/bin/env python
# print __doc__

import sys
import os.path
import threading
import common.adni_tools as adni
import ireg_linear

if len( sys.argv ) < 5:
    print 'Usage: run_ireg_baseline_linear.py <threads> <study> <field_strength> <viscode>'
    exit()
  
nr_threads = int( sys.argv[1] )
study = sys.argv[2]
fs = sys.argv[3]
viscode = sys.argv[4]

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', study )
rreg_params = '/vol/biomedic/users/aschmidt/ADNI/scripts/registration/params-ireg-rigid.txt'
areg_params = '/vol/biomedic/users/aschmidt/ADNI/scripts/registration/params-ireg-affine.txt'

mask_folder = os.path.join( data_folder, 'native/masks_brain' )

output_folder = adni.make_dir( data_folder, 'baseline_linear' )
output_folder_img = adni.make_dir( output_folder, 'images_unstripped' )
output_folder_dof = adni.make_dir( output_folder, 'dof' )

baseline_folder = os.path.join( data_folder, 'native/images_unstripped' )
followup_folder = os.path.join( data_folder, 'native/images_unstripped' )
baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, study, fs, viscode )

class RegistrationThread(threading.Thread):
    def __init__(self, index):
        threading.Thread.__init__(self)
        self.index = index
    def run(self):
        target = baseline_files[self.index]
        source = followup_files[self.index]
        source_base = os.path.basename( source )
        
        out_dof = os.path.join( output_folder_dof, source_base.replace('.nii.gz', '.dof.gz') )
        out_warped = os.path.join( output_folder_img, source_base )

        ireg_linear.run( source, target, 'none', out_dof, rreg_params, areg_params, out_warped )

print 'Found ' + str(len( baseline_files )) + ' image pairs...'
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
    
