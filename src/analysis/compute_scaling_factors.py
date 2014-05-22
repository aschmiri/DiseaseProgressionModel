#! /usr/bin/env python
# print __doc__

import argparse
import csv
import os
from subprocess import check_output
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( '--folder', type=str, default='/vol/medic01/users/aschmidt/projects/Data/ADNI/data/ADNIX/MNI152_linear/dof', help='the path with the dof files with the affine transformations to MNI space' )
parser.add_argument( '--output', type=str, default='./scaling_factors.csv', help='the output file with the scaling factors' )
a = parser.parse_args()

exec_factors = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/affineScaling'

#
# Write a .csv file with the scaling factors
#
def compute_for_study( study ):
    dof_folder = a.folder.replace( 'ADNIX', study )
    files, rids, _, _, _ = adni.read_list_all_data( dof_folder, diagnosis='ALL', study=study, viscode='bl' )
    
    print 'Found', str(len( files )), 'images for', study, '...'
    with open( a.output, 'a' ) as csvfile:
        writer = csv.writer( csvfile, delimiter=',' )
        writer.writerow( ['RID', 'Factor', 'Filename'] )
        for rid, dof in zip( rids, files ):
            factor = float( check_output([exec_factors, dof]) )
            writer.writerow( [rid, factor, dof] )
            #print rid, factor, dof

os.remove( a.output )
compute_for_study( 'ADNI1' )
compute_for_study( 'ADNIGO' )
compute_for_study( 'ADNI2' )
