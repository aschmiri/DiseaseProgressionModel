#! /usr/bin/env python
# print __doc__

import csv
import os.path
import numpy as np
from subprocess import call

rview = '/vol/medic01/users/as12312/Code/Build/merapi/irtk-testing/bin/rview'

datafile_bl  = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/data_bl.csv'
datafile_m24 = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/data_m24.csv'
image_folder = '/vol/biomedic/users/aschmidt/ADNI/data/ADNI1/MNI152_linear/images'

def read_data( datafile, diagnoses ):
    images = []
    virtps = []
    
    with open( datafile, 'rb' ) as csvfile:
        reader = csv.reader( csvfile, delimiter=',' )
        headers = reader.next()
        for row in reader:
            dx = row[headers.index('DX.bl')]
            dof = row[headers.index('FILE')]
            basename = os.path.basename( dof ).replace( '.dof.gz', '.nii.gz' )
            image = os.path.join( image_folder, basename )
            if dx in diagnoses and os.path.exists( image ):
                images.append( image )
                virtps.append( float( row[headers.index('VIRT_P')] ) )
            
    return (images, virtps)


images, virtp = read_data( datafile_bl, ['AD'] )
images = np.array( images )
virtp = np.array( virtp )

indices = np.argsort( virtp )
indices = indices[::-1]

images_sorted = images[ indices ]
virtp_sorted = virtp[ indices ]

for i in range( len( images_sorted ) ):
    print 'VIRT_P: ' + str(virtp_sorted[i]) + ' ' + images_sorted[i]
    call([rview, images_sorted[i]])

