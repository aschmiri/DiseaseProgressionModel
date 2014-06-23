#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import threading
from subprocess import call
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( '-n', '--nr_threads', dest = 'nr_threads', type=int, default = 1 )
parser.add_argument( '-s', '--spacing', dest = 'sx', type=str, default = '10' )
a = parser.parse_args()

data_folder = os.path.join( adni.data_folder, a.study )
mask_folder = os.path.join( data_folder, 'native/masks_brain' )

dof_folder_lin    = os.path.join( data_folder, 'baseline_linear', 'dof' )
dof_folder_nonlin = os.path.join( data_folder, 'baseline_' + a.trans + '_' + a.sx + 'mm_after_linear', 'dof' )

seg_folder     = os.path.join( '/vol/medic01/users/aschmidt/projects/Data/ADNI/data', a.study, 'native' )
seg_folder_in  = os.path.join( seg_folder, 'seg_138regions_baseline' )
seg_folder_out = adni.make_dir( seg_folder, 'seg_138regions_followup_' + a.trans + '_' + a.sx + 'mm' )

baseline_folder = os.path.join( data_folder, 'native/images_unstripped' )
followup_folder = os.path.join( data_folder, 'native/images_unstripped' )
baseline_files, followup_files = adni.get_baseline_and_followup( baseline_folder, followup_folder, a.study, a.viscode )

class RegistrationThread(threading.Thread):
    def __init__(self, index):
        threading.Thread.__init__(self)
        self.index = index
    def run(self):
        target = baseline_files[self.index]
        target_base = os.path.basename( target )
        study_bl = adni.detect_study( target )
        source = followup_files[self.index]
        source_base = os.path.basename( source )
        
        dof_lin    = os.path.join( dof_folder_lin, source_base.replace('.nii.gz', '.dof.gz') )
        dof_nonlin = os.path.join( dof_folder_nonlin, source_base.replace('.nii.gz', '.dof.gz') )
         
        target_seg = os.path.join( seg_folder_in.replace(a.study, study_bl), 'EM-' + target_base )
        out_seg    = os.path.join( seg_folder_out, source_base )
         
        if not os.path.isfile( target_seg ):
            print 'Segmentation ' + target_seg + ' does not exists!'
        if not os.path.isfile( dof_nonlin ):
            print 'DOF file ' + dof_nonlin + ' does not exists!'
        elif os.path.isfile( out_seg ):
            print 'File ' + out_seg + ' already exists!'
        else:
            #
            # Run transform masks
            #
            print '--------------------'
            print 'Starting: transformation'
            print 'Target:     ' + target
            print 'Source:     ' + source
            print 'DOF lin:    ' + dof_lin
            print 'DOF nonlin: ' + dof_nonlin
            print 'Seg in:     ' + target_seg
            print 'Seg out:    ' + out_seg
            
            dof_lin_ffd = out_seg.replace( '.nii.gz', '_lin.dof.gz' )
            dof_combined = out_seg.replace( '.nii.gz', '_combined.dof.gz' )
            
            call([ 'ffdcreate', dof_lin_ffd, '-dofin', dof_lin ])
            call([ 'ffdcompose', dof_nonlin, dof_lin_ffd, dof_combined ])
            call([ 'transformation', target_seg, out_seg, '-dofin', dof_combined, '-target', target, '-nn', '-matchInputType', '-invert' ])
            
#             call([ 'transformation', target_seg, out_seg_nonlin, '-dofin', dof_nonlin, '-target', target, '-nn', '-matchInputType', '-invert' ])
#             call([ 'transformation', out_seg_nonlin, out_seg, '-dofin', dof_lin, '-target', source, '-nn', '-matchInputType', '-invert' ])
            
            call([ 'rm', dof_lin_ffd ])
            call([ 'rm', dof_combined ])
            
            print '--------------------'
        

print 'Found', str(len( baseline_files )), 'image pairs...'
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
