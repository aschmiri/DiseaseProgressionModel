#! /usr/bin/env python
# print __doc__

import os.path
import argparse
from subprocess import call
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...' )
a = parser.parse_args()

execAverage = 'ffdaverage'

base_folder = '/vol/biomedic/users/aschmidt/ADNI'
dof_folder_adni1 = os.path.join( base_folder, 'data/ADNI1/MNI152_svffd_10mm_followup_to_baseline/dof' )
dof_folder_adni2 = os.path.join( base_folder, 'data/ADNI2/MNI152_svffd_10mm_followup_to_baseline/dof' )

velocities = adni.get_baseline_transformations( dof_folder_adni1, dof_folder_adni2, a.viscode, a.diagnosis )
out_average_velo = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_0/velo_' + a.viscode + '_' + a.diagnosis + '.dof.gz'    

print 'Found ' + str(len( velocities )) + ' velocities in total for viscode ' + a.viscode + '...'

if os.path.isfile( out_average_velo ):
    print 'Image already exists: ' + out_average_velo
else:
    print '--------------------'
    print 'Output velocities: ' + out_average_velo
      
    call([ execAverage, out_average_velo] + velocities )
