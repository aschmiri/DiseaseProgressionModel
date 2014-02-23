#! /usr/bin/env python
# print __doc__

import csv
import matplotlib.pyplot as plt
import numpy as np

datafile_bl  = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_1/data_m24_AD_fielddist.csv'

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
states = [states_cn, states_mci , states_ad]

print len(states_ad)

num_bins = 41

y, bin_edges = np.histogram( states_ad, num_bins )
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
#n, bins, patches = plt.hist(states, num_bins, normed=1, color=['green','yellow','red'], alpha = 0.5)
n, bins, patches = plt.hist(states_ad, num_bins, normed=1, color='red', alpha = 0.5, histtype='step' )
#plt.plot( bin_centers, y, '-' )
plt.xlabel('Virtual disease state')
plt.ylabel('Probability')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()