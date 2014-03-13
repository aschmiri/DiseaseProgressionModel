#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import threading
from subprocess import call
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'field_strength', type=str,  help='the field strength, usually 1.5 for ADNI1 and 3 otherwise' )
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( '-n', '--nr_threads', dest = 'nr_threads', type=int, default = 1 )
parser.add_argument( '-s', '--spacing', dest = 'sx', type=str, default = '10' )
a = parser.parse_args()

target_mni = os.path.join( adni.mni_folder, 'MNI152_T1_1mm_brain.nii' )

data_folder = os.path.join( adni.data_folder, a.study )
baseline_folder = os.path.join( data_folder, 'native/images' )
followup_folder = os.path.join( data_folder, 'baseline_linear/images' )

mni_linear_folder = os.path.join( data_folder, 'MNI152_linear' )
mni_linear_folder_dof = os.path.join( mni_linear_folder, 'dof' )

mni_nonlin_folder = os.path.join( data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear' )
mni_nonlin_folder_dof = os.path.join( mni_nonlin_folder, 'dof' )

out_folder_linear = adni.make_dir( data_folder, 'MNI152_linear_via_baseline' )
out_folder_linear_img = adni.make_dir( out_folder_linear, 'images' )
out_folder_nonlin_img = adni.make_dir( mni_nonlin_folder, 'images_' + a.viscode )

baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, a.study, a.field_strength, a.viscode )


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
        if a.study == 'ADNI2' and baseline.find( 'ADNIGO' ) > -1:
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
     
    if thread_ctr == a.nr_threads:
        for t in threads:
            t.join()
        threads = []
        thread_ctr = 0