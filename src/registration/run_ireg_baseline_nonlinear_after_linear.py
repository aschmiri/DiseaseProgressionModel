#! /usr/bin/env python
# print __doc__

import sys
import os.path
import threading
import common.adni_tools as adni
import ireg_nonlinear

if len( sys.argv ) < 7:
    print 'Usage: run_ireg_baseline_nonlinear_after_linear.py <threads> <study> <field_strength> <viscode> <transformation> <spacing>'
    exit()
  
nr_threads = int( sys.argv[1] )
study = sys.argv[2]
fs = sys.argv[3]
viscode = sys.argv[4]
trans = sys.argv[5]
sx = sys.argv[6]

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', study )

ireg_params = '/vol/biomedic/users/aschmidt/ADNI/scripts/registration/params-ireg-' + trans + '-' + sx + 'mm.txt'

baseline_folder = os.path.join( data_folder, 'native/images' )
followup_folder = os.path.join( data_folder, 'baseline_linear/images' )
baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, study, fs, viscode )

output_folder = adni.make_dir( data_folder, 'baseline_' + trans + '_' + sx + 'mm_after_linear' )
output_folder_img = adni.make_dir( output_folder, 'images' )
output_folder_dof = adni.make_dir( output_folder, 'dof' )

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
        
        ireg_nonlinear.run( source, target, 'none', out_dof, ireg_params, out_warped )

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
