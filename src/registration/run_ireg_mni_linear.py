#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import threading
import common.adni_tools as adni
import ireg_linear

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'field_strength', type=str,  help='the field strength, usually 1.5 for ADNI1 and 3 otherwise' )
parser.add_argument( '-n', '--nr_threads', dest = 'nr_threads', type=int, default = 1 )
a = parser.parse_args()

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', a.study )

target_mni = '/vol/medic01/users/aschmidt/projects/Data/MNI152-Template/MNI152_T1_1mm_brain.nii'
rreg_params = '/vol/biomedic/users/aschmidt/ADNI/scripts/registration/params-ireg-rigid.txt'
areg_params = '/vol/biomedic/users/aschmidt/ADNI/scripts/registration/params-ireg-affine.txt'

baseline_folder = os.path.join( data_folder, 'native/images' )
baseline_files = adni.get_baseline( baseline_folder, a.study, a.field_strength )

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
     
    if thread_ctr == a.nr_threads:
        for t in threads:
            t.join()
        threads = []
        thread_ctr = 0
