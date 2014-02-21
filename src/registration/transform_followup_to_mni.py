#! /usr/bin/env python
# print __doc__

import sys
import os.path
import threading
from subprocess import call
import common.adni_tools as adni

if len( sys.argv ) < 7:
    print 'Usage: tarnsform_baseline_to_mni.py <threads> <study> <field_strength> <viscode> <transformation> <spacing>'
    exit()
  
nr_threads = int( sys.argv[1] )
study = sys.argv[2]
fs = sys.argv[3]
viscode = sys.argv[4]
trans = sys.argv[5]
sx = sys.argv[6]

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
data_folder = os.path.join( base_folder, 'data', study )

target_mni = '/vol/medic01/users/aschmidt/projects/Data/MNI152-Template/MNI152_T1_1mm_brain.nii'

baseline_folder = os.path.join( data_folder, 'native/images' )
followup_folder = os.path.join( data_folder, 'baseline_linear/images' )

mni_linear_folder = os.path.join( data_folder, 'MNI152_linear' )
mni_linear_folder_dof = os.path.join( mni_linear_folder, 'dof' )

mni_nonlin_folder = os.path.join( data_folder, 'MNI152_' + trans + '_' + sx + 'mm_after_linear' )
mni_nonlin_folder_dof = os.path.join( mni_nonlin_folder, 'dof' )

out_folder_linear = adni.make_dir( data_folder, 'MNI152_linear_via_baseline' )
out_folder_linear_img = adni.make_dir( out_folder_linear, 'images' )
out_folder_nonlin_img = adni.make_dir( mni_nonlin_folder, 'images_' + viscode )

baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, study, fs, viscode )


class RegistrationThread(threading.Thread):
    def __init__(self, index):
        threading.Thread.__init__(self)
        self.index = index
    def run(self):
        baseline = baseline_files[self.index]
        followup = followup_files[self.index]
        baseline_base = os.path.basename( baseline )
        followup_base = os.path.basename( followup )
        aff = os.path.join( mni_linear_folder_dof, baseline_base.replace('.nii.gz', '.dof.gz') )
        dof = os.path.join( mni_nonlin_folder_dof, baseline_base.replace('.nii.gz', '.dof.gz') )
        
        # Rename files if baseline scan is found in ADNIGO cohort.
        if study == 'ADNI2' and baseline.find( 'ADNIGO' ) > -1:
            aff = aff.replace( 'ADNI2', 'ADNIGO' )
            dof = dof.replace( 'ADNI2', 'ADNIGO' )

        out_image_affine = os.path.join( out_folder_linear_img, followup_base )
        out_image_nonlin = os.path.join( out_folder_nonlin_img, followup_base )
        
        if os.path.isfile( out_image_affine ):
            print 'Image already exists: ' + out_image_nonlin
        else:
            print '--------------------'
            print 'In image:  ' + followup
            print 'Affine:    ' + aff
            print 'Out image: ' + out_image_affine
            call([ 'transformation', followup, out_image_affine, '-dofin', aff, '-target', target_mni, '-cspline', '-matchInputType', '-Sp', '0' ])
        
        if os.path.isfile( out_image_nonlin ):
            print 'Image already exists: ' + out_image_nonlin
        else:     
            print '--------------------'
            print 'In image:  ' + out_image_affine
            print 'Nonlinear: ' + dof
            print 'Out image: ' + out_image_nonlin
            call([ 'transformation', out_image_affine, out_image_nonlin, '-dofin', dof, '-target', target_mni, '-cspline', '-matchInputType', '-Sp', '0' ])
                 
print 'Found ' + str(len( followup_files )) + ' images...'
thread_ctr = 0
threads = []
for i in range( len( followup_files ) ):
    thread = RegistrationThread(i)
    thread.start()
    threads.append(thread)
    thread_ctr += 1
     
    if thread_ctr == nr_threads:
        for t in threads:
            t.join()
        threads = []
        thread_ctr = 0