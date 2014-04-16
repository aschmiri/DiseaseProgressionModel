#! /usr/bin/env python
# print __doc__

import csv
import matplotlib.pyplot as plt
import numpy as np

datafile_m12  = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_0_longit/data_m12_AD.csv'
datafile_m24  = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_0_longit/data_m24_AD.csv'

def read_data( datafile, diagnoses ):
    states = []
    
    with open( datafile, 'rb' ) as csvfile:
        reader = csv.reader( csvfile, delimiter=',' )
        headers = reader.next()
        for row in reader:
            dx = row[headers.index('DX.bl')]
            if dx in diagnoses:
                states.append( float( row[headers.index('DPI')] ) )
            
    return np.array( states )

#
# Read data
states_m12_ad  = read_data( datafile_m12, ['AD'] )
states_m24_ad  = 2 * read_data( datafile_m24, ['AD'] )

plt.scatter( states_m24_ad, states_m12_ad, color=(0,0.38,0.48), alpha=0.2 )

sxx=0
syy=0
sxy=0
for i in range(len(states_m24_ad)):
    sx = states_m24_ad[i]
    sy = states_m12_ad[i]
    sxx += sx * sx
    syy += sy * sy
    sxy += sx * sy


f = sxy/sxx
print f
x = np.arange(0, 22)
y = x * f

e=0
for i in range(len(states_m24_ad)):
    e += np.square( states_m12_ad[i] - f * states_m24_ad[i] )
e = np.sqrt(e/len(states_m24_ad))
print e
    
    
plt.plot( x, y, color=(0,0.38,0.48) )


plt.plot( [0,21], [0,21], '--', color='grey' )

plt.xlabel('Virtual disease state, m24 model')
plt.ylabel('Virtual disease state, m12 model')
plt.xlim([-1,21])
plt.ylim([-1,21])

# Tweak spacing to prevent clipping of ylabel
#plt.subplots_adjust(left=0.15)
plt.show()