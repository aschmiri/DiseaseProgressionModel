#! /usr/bin/env python
# print __doc__

import numpy as np
import argparse
import csv
import os.path
from subprocess import call
import common.atlas_tools as at

parser = argparse.ArgumentParser()
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...' )
parser.add_argument( '-i', '--iteration', dest='iteration', type=int, default=1 )
parser.add_argument( '-r', '--required_subjects', dest='required_subjects', type=int, default=50 )
parser.add_argument( '--min', dest='state_min', type=float, default=0 )
parser.add_argument( '--max', dest='state_max', type=float, default=10 )
parser.add_argument( '--steps', dest = 'state_steps', type=int, default=11 )
a = parser.parse_args()

exec_atlas = 'atlas'
exec_ffdaverage = 'ffdaverage'
exec_transform = 'transformation'

method_folder = 'MNI152_svffd_10mm_after_linear'
type_folder_images = 'images_normalised'
type_folder_dofs = 'dof'

atlas_folder = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_' + str(a.iteration)
datafile = os.path.join( atlas_folder, 'data_' + a.viscode + '_' + a.diagnosis + '.csv' )

_, _, _, states, images = at.read_datafile( datafile, a.diagnosis )

for state in np.linspace( a.state_min, a.state_max, a.state_steps ):
    # Find sigma and corresponding images
    sigma, weights, indices = at.adaptive_kernel_regression( states, state, required_subjects = a.required_subjects )
    
    selected_images = images[indices]
    selected_weights = weights[indices]
    
    # Define the base name of output files
    atlas_base = 'atlas_' + a.viscode + '_' + a.diagnosis + '_' + str(state)
    
    # Print data file for IRTK ffd averaging
    data_file_dofs = os.path.join( atlas_folder, atlas_base + '_dofs.txt' )
    with open( data_file_dofs, 'wb') as csvfile:
        csv_writer = csv.writer( csvfile, delimiter=' ' )
        for i in range(len( selected_images )):
            dof = at.find_file( selected_images[i], method_folder, type_folder_dofs )
            csv_writer.writerow([ dof, selected_weights[i] ])
     
    # Print data file for IRTK image averaging
    data_file_images = os.path.join( atlas_folder, atlas_base + '_images.txt' )
    with open( data_file_images, 'wb') as csvfile:
        csv_writer = csv.writer( csvfile, delimiter=' ' )
        for i in range(len( selected_images )):
            image = at.find_file( selected_images[i], method_folder, type_folder_images )
            csv_writer.writerow([ image, selected_weights[i] ])           
            
    
    # Call IRTK 'ffdaverage'
    out_average_dof = os.path.join( atlas_folder, atlas_base + '_average_ffd.dof.gz' )
    if not os.path.exists( out_average_dof ):
        print '--------------------'
        print 'Starting: ffdaverage'
        print 'Data:   ' + data_file_dofs
        print 'Output: ' + out_average_dof
        call( [exec_ffdaverage, out_average_dof, '-dofnames', data_file_dofs ] )
    
    # Call IRTK 'atlas'
    out_average_image = os.path.join( atlas_folder, atlas_base + '_average_image.nii.gz' )
    if not os.path.exists( out_average_image ):
        print '--------------------'
        print 'Starting: ffdaverage'
        print 'Data:   ' + data_file_images
        print 'Output: ' + out_average_image
        call( [exec_atlas, out_average_image, '-imagenames', data_file_images ] )
    
    # Transform average image
    out_average_image_transformed = os.path.join( atlas_folder, atlas_base + '_average_image_transformed.nii.gz' )
    if not os.path.exists( out_average_image_transformed ):
        print '--------------------'
        print 'Starting: transformation'
        print 'Input:  ' + out_average_image
        print 'DOF:    ' + out_average_dof
        print 'Output: ' + out_average_image_transformed
        call([ exec_transform, out_average_image, out_average_image_transformed, 
               '-dofin', out_average_dof, '-invert', '-cspline' ])

