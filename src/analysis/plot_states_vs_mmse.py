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
                states.append( float( row[headers.index('VIRT_NMI')] ) )
            
    return (ages, mmses, states)

#
# Read data
ages_cn,  mmses_cn,  states_cn  = read_data( datafile_bl, ['CN'] )
ages_mci, mmses_mci, states_mci = read_data( datafile_bl, ['EMCI','LMCI'] )
ages_ad,  mmses_ad,  states_ad  = read_data( datafile_bl, ['AD'] )

print 'Pearson CC: Disease state // MMSE'
print 'CN:  ', pearsonr(mmses_cn, states_cn)
print 'MCI: ', pearsonr(mmses_mci, states_mci)
print 'AD:  ', pearsonr(mmses_ad, states_ad)

print 'Pearson CC: Age // MMSE'
print 'CN:  ', pearsonr(ages_cn, states_cn)
print 'MCI: ', pearsonr(ages_mci, states_mci)
print 'AD:  ', pearsonr(ages_ad, states_ad)

f, (ax1,ax2) = plt.subplots(1,2)
#ax1.scatter( states_cn,  mmses_cn,  color='green' )
#ax1.scatter( states_mci, mmses_mci, color='yellow' )
ax1.scatter( states_ad,  mmses_ad,  color='red' )

#ax2.scatter( states_cn,  ages_cn,  color='green' )
#ax2.scatter( states_mci, ages_mci, color='yellow' )
ax2.scatter( states_ad,  ages_ad,  color='red' )

plt.show()
