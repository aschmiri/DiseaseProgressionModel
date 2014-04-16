#! /usr/bin/env python
# print __doc__

import csv
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

datafile_bl  = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_0/data_sym_m24_AD.csv'

def read_data( datafile, diagnoses ):
    ages = []
    mmses = []
    states = []
    
    with open( datafile, 'rb' ) as csvfile:
        reader = csv.reader( csvfile, delimiter=',' )
        headers = reader.next()
        for row in reader:
            dx = row[headers.index('DX.bl')]
            if dx in diagnoses:
                ages.append( float( row[headers.index('AGE')] ) )
                mmses.append( int( row[headers.index('MMSE')] ) )
                states.append( float( row[headers.index('DPI')] ) )
            
    return (ages, mmses, states)

#
# Read data
ages_cn,  mmses_cn,  states_cn  = read_data( datafile_bl, ['CN'] )
ages_mci, mmses_mci, states_mci = read_data( datafile_bl, ['EMCI','LMCI'] )
ages_ad,  mmses_ad,  states_ad  = read_data( datafile_bl, ['AD'] )

print 'Pearson CC: Disease state // MMSE'
#print 'CN:  ', pearsonr(mmses_cn, states_cn)
#print 'MCI: ', pearsonr(mmses_mci, states_mci)
print 'AD:  ', pearsonr(mmses_ad, states_ad)

print 'Pearson CC: Disease state // Age'
#print 'CN:  ', pearsonr(ages_cn, states_cn)
#print 'MCI: ', pearsonr(ages_mci, states_mci)
print 'AD:  ', pearsonr(ages_ad, states_ad)

f, (ax1,ax2) = plt.subplots(1,2)
#ax1.scatter( states_cn,  mmses_cn,  color='green' )
#ax1.scatter( states_mci, mmses_mci, color='yellow' )
#ax1.scatter( states_ad,  mmses_ad,  color='red' )

#ax2.scatter( states_cn,  ages_cn,  color='green' )
#ax2.scatter( states_mci, ages_mci, color='yellow' )
#ax2.scatter( states_ad,  ages_ad,  color='red' )


import numpy as np
import scipy.stats as stats
states_ad = np.array( states_ad )
ages_ad = np.array( ages_ad )
mean_state = np.mean( states_ad )
slope, intercept, _, _, _ = stats.linregress( ages_ad , states_ad )
states_ad_regressed  = states_ad - (ages_ad * slope + intercept) + mean_state
print 'Pearson CC: Disease state // MMSE'
print 'AD:  ', pearsonr(mmses_ad, states_ad_regressed)

print 'Pearson CC: Disease state // Age'
print 'AD:  ', pearsonr(ages_ad, states_ad_regressed)

ax1.scatter( ages_ad, states_ad,  color='red' )
ax2.scatter( ages_ad, states_ad_regressed,  color='red' )
x = [2,6]
ax1.plot((x - intercept) / slope,x,'-b' )
ax1.set_ylim([-2,16])
ax1.set_xlim([50,95])
ax2.set_ylim([-2,16])
ax2.set_xlim([50,95])

plt.show()
