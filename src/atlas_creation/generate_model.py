#! /usr/bin/env python
# print __doc__

import os.path
import argparse
from subprocess import call

parser = argparse.ArgumentParser( formatter_class=argparse.ArgumentDefaultsHelpFormatter )
parser.add_argument( 'viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...' )
parser.add_argument( 'diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...' )
parser.add_argument( '-i', '--iteration', dest='iteration', type=int, default=0, help='the current iteration' )
parser.add_argument( '-r', '--required_subjects', dest='required_subjects', type=int, default=50, help='number of required subjects' )
parser.add_argument( '--state_min', type=float, default=0, help='minimal state' )
parser.add_argument( '--state_max', type=float, default=15, help='maximal state' )
parser.add_argument( '--state_stepsize', type=float, default = 0.25, help='step size for model generation' )
a = parser.parse_args()

exec_model = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/model_generation'

base_folder = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_' + str(a.iteration)
mean_atlas = '/vol/medic01/users/aschmidt/projects/Data/MNI152-Template/MNI152_T1_1mm_brain2.nii'
model_prefix = os.path.join( base_folder, 'model_' + a.viscode + '_' + a.diagnosis + '_s' )

linear = True if a.iteration == 0 else False
if linear:
    velo = os.path.join( base_folder, 'velo_' + a.viscode + '_' + a.diagnosis + '.dof.gz' )
    call([ exec_model, mean_atlas, velo, 
          '-min', str(a.state_min), '-max', str(a.state_max), '-stepsize', str(a.state_stepsize),
          '-saveimg', '-prefix', model_prefix, '-linear'])
else:
    velo_file = os.path.join( base_folder, 'velo_' + a.viscode + '_' + a.diagnosis + '.txt' )
    call([ exec_model, mean_atlas, velo_file, 
          '-min', str(a.state_min), '-max', str(a.state_max), '-stepsize', str(a.state_stepsize),
          '-saveimg', '-prefix', model_prefix ])
    