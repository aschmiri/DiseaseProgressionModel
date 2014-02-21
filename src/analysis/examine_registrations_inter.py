#! /usr/bin/env python
# print __doc__

import argparse
import os.path
from subprocess import call
import common.atlas_tools as at

parser = argparse.ArgumentParser()
parser.add_argument( 'trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic' )
parser.add_argument( 'state', type=float, help='the state for which relevant images should be registered' )
parser.add_argument( '-s', '--spacing', dest='sx', type=str, default='10' )
parser.add_argument( '-r', '--required_subjects', dest='required_subjects', type=int, default=20 )
parser.add_argument( '-i', '--rid', dest='rid', type=str, default=None )
a = parser.parse_args()
    
rview = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/rview'

datafile = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_1/data_m24_AD.csv'
rids, _, _, states, images = at.read_all_data( datafile, ['AD'] )

sigma, weights, indices = at.adaptive_kernel_regression( states, a.state, required_subjects = a.required_subjects )

selected_rids = rids[indices]
selected_images = images[indices]
selected_weights = weights[indices]

data_folder = '/vol/biomedic/users/aschmidt/ADNI/data/ADNI'
dof_folder = os.path.join( data_folder, 'MNI152_intra_' + a.trans + '_' + a.sx + 'mm', 'dof' )

print 'Found ' + str(len( selected_images )) + ' relevant images for state ' + str(a.state) + '...'
for i in range( len( selected_images ) ):
    for j in range( len( selected_images ) ):
        target = selected_images[i]
        source = selected_images[j]
        target_rid = selected_rids[i]
        source_rid = selected_rids[j]
        dof_basename = str(source_rid) + '_to_' + str(target_rid)

        dof = os.path.join( dof_folder, dof_basename + '.dof.gz' )
    
        if os.path.exists( dof ):
            if a.rid == None or a.rid == target_rid:
                print '--------------------'
                print 'Target: ' + target
                print 'Source: ' + source
                print 'DOF:    ' + dof
                
                call([ rview, target, source, dof, '-res', '1.5', '-mix' ])
