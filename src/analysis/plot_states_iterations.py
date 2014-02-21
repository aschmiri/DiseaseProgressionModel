#! /usr/bin/env python
# print __doc__

import csv
import matplotlib.pyplot as plt
import numpy as np

datafile_m0  = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_0/data_m24_AD.csv'
datafile_m1  = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_1/data_m24_AD.csv'

def read_data( datafile, diagnoses ):
    states = []
    
    with open( datafile, 'rb' ) as csvfile:
        reader = csv.reader( csvfile, delimiter=',' )
        headers = reader.next()
        for row in reader:
            dx = row[headers.index('DX.bl')]
            if dx in diagnoses:
                states.append( float( row[headers.index('VIRT_NMI')] ) )
            
    return states

#
# Read data
states_m0_cn  = read_data( datafile_m0, ['CN'] )
states_m0_mci  = read_data( datafile_m0, ['EMCI', 'LMCI'] )
states_m0_ad  = read_data( datafile_m0, ['AD'] )

states_m1_cn  = read_data( datafile_m1, ['CN'] )
states_m1_mci  = read_data( datafile_m1, ['EMCI', 'LMCI'] )
states_m1_ad  = read_data( datafile_m1, ['AD'] )

plt.scatter( states_m0_ad, states_m1_ad, color='yellow' )

coefs =  np.polyfit( states_m0_ad, states_m1_ad, 2)
print coefs
x = np.arange(0, 10)
y = np.square(x) * coefs[0] + x * coefs[1] + coefs[2]
#y = x * coefs[0] + coefs[1]
plt.plot( x, y )


plt.xlabel('Virtual disease state m0')
plt.ylabel('Virtual disease state m1')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()