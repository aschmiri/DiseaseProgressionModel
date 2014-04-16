#! /usr/bin/env python
# print __doc__

import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

datafile_bl  = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/data_mni_bl.csv'
datafile_m24 = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/data_mni_m12.csv'

def read_data( datafile ):
    rids = []
    states = []
    
    with open( datafile, 'rb' ) as csvfile:
        reader = csv.reader( csvfile, delimiter=',' )
        headers = reader.next()
        for row in reader:
            dx = row[headers.index('DX.bl')]
            if dx == 'AD':
                rids.append( int( row[headers.index('RID')] ) )
                states.append( float( row[headers.index('DPI')] ) )
            
    return (rids, states)

#
# Read data
rids_bl,  states_bl_unsorted  = read_data( datafile_bl )
rids_m24, states_m24_unsorted = read_data( datafile_m24 )

#
# Sort using the RIDs
states_bl  = []
states_m24 = []
for i in range( len(rids_bl) ):
    baseline_rid = rids_bl[i]
    for j in range( len(rids_m24) ):
        followup_rid = rids_m24[j]
        if baseline_rid == followup_rid:
            states_bl.append( states_bl_unsorted[i] )
            states_m24.append( states_m24_unsorted[j] )

#
# Calculatettest
t, prob = scipy.stats.ttest_rel( states_bl, states_m24 )
print 'T    ' + str(t)
print 'P    ' + str(prob)

#
# Calculate differences
states_diff = []
for i in range( len(states_bl) ):
    diff = states_m24[i] - states_bl[i]
    states_diff.append( diff )
    
print 'Min  ' + str(np.amin(states_diff))
print 'Max  ' + str(np.amax(states_diff))
print 'Mean ' + str(np.mean(states_diff))

#
# Plot
for i in range(len(states_bl)):
    plt.plot((0,1),(states_bl[i],states_m24[i]) , color='red')
plt.plot( (0,1), (np.mean(states_bl),np.mean(states_m24)), color = 'black', lw = 3 )

plt.show()
