#! /usr/bin/env python
# print __doc__

import csv
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import common.atlas_tools as at

datafile_bl  = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_1/data_m24_AD.csv'

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
states_cn  = read_data( datafile_bl, ['CN'] )
states_mci  = read_data( datafile_bl, ['EMCI', 'LMCI'] )
states_ad  = read_data( datafile_bl, ['AD'] )

x = np.arange(0, 10.5, 0.5)
y = np.zeros(21)
for i in range(len(states_ad)):
    a = np.floor( states_ad[i] )
    b = states_ad[i] - a
    if b == 0:
        index = 2 * a
        y[index] += 1
    elif b == 0.5:
        index = 2 * a + 1
        y[index] += 1
    elif b == 0.25:
        index1 = 2 * a
        index2 = 2 * a + 1 
        y[index1] += 0.5
        y[index2] += 0.5
    elif b == 0.75:
        index1 = 2 * a + 1
        index2 = 2 * a + 2 
        y[index1] += 0.5
        y[index2] += 0.5
    else:
        print ' unknown value '

plt.bar( x, y, align='center', width=0.4, linewidth=0, color=(0,0.38,0.48), alpha=0.7 )
plt.xlabel('Virtual disease state')
plt.ylabel('Probability')
plt.xlim([-0.5,10.5])
# ocean green (0,0.29,0.35) // '#004B5A'
prob_x = np.arange(-1, 11, 0.05)
for i in range(11):
    sigma, _, _ = at.adaptive_kernel_regression( np.array( states_ad ), i, required_subjects=50 )
    prob_y = mlab.normpdf( prob_x, i, sigma )
    plt.plot( prob_x, 6*prob_y, linewidth=1.5, color='black' )
    
    
# Tweak spacing to prevent clipping of ylabel
#plt.subplots_adjust(left=0.15)
plt.show()
    









    