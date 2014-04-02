#! /usr/bin/env python
# print __doc__

import os.path
import argparse
from subprocess import call
import numpy as np
import csv
import common.adni_tools as adni
import common.atlas_tools as at

parser = argparse.ArgumentParser()
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...' )
parser.add_argument( '-i', '--iteration', dest='iteration', type=int, default=1 )
parser.add_argument( '-r', '--required_subjects', dest='required_subjects', type=int, default=20 )
parser.add_argument( '-t', '--trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( '--min', dest='state_min', type=float, default=0 )
parser.add_argument( '--max', dest='state_max', type=float, default=10 )
parser.add_argument( '--steps', dest='state_steps', type=int, default=11 )
parser.add_argument( '--postfix', type=str, default='' )
a = parser.parse_args()

execAverage = 'ffdaverage'

dof_folder_adni1 = os.path.join( adni.data_folder, 'ADNI1', 'MNI152_' + a.trans +'_10mm_followup_to_baseline/dof' )
dof_folder_adni2 = os.path.join( adni.data_folder, 'ADNI2', 'MNI152_' + a.trans +'_10mm_followup_to_baseline/dof' )

velocities = adni.get_baseline_transformations( dof_folder_adni1, dof_folder_adni2, a.viscode, a.diagnosis )

print 'Found ' + str(len( velocities )) + ' velocities in total for viscode ' + a.viscode + '...'

atlas_folder = os.path.join( adni.project_folder, 'atlas' )
atlas_folder_in = os.path.join( atlas_folder, 'model_' + str(a.iteration-1) + a.postfix )
atlas_folder_out = os.path.join( atlas_folder, 'model_' + str(a.iteration) + a.postfix )
input_datafile =  os.path.join( atlas_folder_in, 'data_' + a.trans +'_' + a.viscode + '_' + a.diagnosis + '.csv' )

#-------------------------------------------------------------------------------
# get data
rids_bl, _, _, states_bl, _ = at.read_datafile( input_datafile, a.diagnosis )
velocities_fu, rids_fu = adni.get_baseline_transformations_and_rids( dof_folder_adni1, dof_folder_adni2, a.viscode, a.diagnosis )

#-------------------------------------------------------------------------------
# sort data to match velocities with templates
dofs = []
states = []

for i in range( len(rids_bl) ):
    baseline_rid = rids_bl[i]
    for j in range( len(rids_fu) ):
        followup_rid = rids_fu[j]
        if baseline_rid == followup_rid:
            states.append( states_bl[i] )
            dofs.append( velocities_fu[j] )
                
dofs = np.array( dofs )
states = np.array( states )
print 'Found ' + str(len( dofs )) + ' velocities with corresponding disease state...'

#-------------------------------------------------------------------------------
# sort data
for state in np.linspace( a.state_min, a.state_max, a.state_steps ):
    # Find sigma and corresponding images
    sigma, weights, indices = at.adaptive_kernel_regression( states, state, required_subjects = a.required_subjects )
    
    selected_dofs = dofs[indices]
    selected_weights = weights[indices]
    selected_states = states[indices]
    
    statestr = '%.2f' % state
    out_average_velo = os.path.join( atlas_folder_out, 'velo_' + a.viscode + '_' + a.diagnosis + '_s' + statestr + ".dof.gz" )
    
    if os.path.isfile( out_average_velo ):
        print 'Image already exists: ' + out_average_velo
    else:
        dof_names_file = os.path.join( atlas_folder_out, 'velo_' + a.viscode + '_' + a.diagnosis + '_s' + statestr + ".txt" )
        with open( dof_names_file, 'wb') as csvfile:
            csv_writer = csv.writer( csvfile, delimiter=' ' )
            for i in range(len( selected_dofs )):
                csv_writer.writerow([ selected_dofs[i], str(selected_states[i]) ])

        print '--------------------'
        print 'DOF names file:    ' + dof_names_file
        print 'Output velocities: ' + out_average_velo
          
        call([ execAverage, out_average_velo, '-dofnames', dof_names_file, '-gaussian', str(state), str(sigma)])

# Save model file required for 'model_generation'
model_file = os.path.join( atlas_folder_out, 'velo_' + a.trans +'_' + a.viscode + '_' + a.diagnosis + '.txt' )
if not os.path.isfile( model_file ):
    call([ 'ls -rt -d -1 ' + atlas_folder_out + '/*' + a.viscode + '_' + a.diagnosis + '*.dof.gz > ' + model_file ], shell=True )
    
    
    
    
