#! /usr/bin/env python
# print __doc__

import os.path
import csv
import numpy as np
import re
from scipy.optimize import curve_fit
import common.adni_tools as adni
import matplotlib.pyplot as plt

def sigmoid( t, t0, r, lower, upper ):
    y = lower + (upper - lower) * ( 1 / (1 + np.exp(-r*(t-t0))) )
    return y

list_r = []
list_t0 = []
list_up = []

def plot_trajectories( traj_x, traj_y, traj_dx ):
    if len( traj_x ) > 5 and traj_dx[-1] == 'AD': # and traj_dx[0] != 'AD':
        # Convert data
        xdata = np.array(traj_x)
        ydata = np.array(traj_y)
        dx = np.array(traj_dx)
        
        try:
            # Fit curve
            popt, _ = curve_fit( sigmoid, xdata, ydata, p0 = (25, 1, min(ydata), max(ydata)) )
            ydata = ydata/popt[2]
            
            if popt[0] > 0 and popt[3] / popt[2] < 2:
                print popt

                # Plot curve
                x = np.linspace( xdata[0], xdata[-1], 100)
                x = np.linspace( 0, 48, 100)
                y = sigmoid( x, *popt) / popt[2]
                plt.plot( x, y, label='fit1', color='b')
                
                list_r.append( np.sign( popt[1] ) * np.log( np.abs ( popt[1] ) ) )
                list_t0.append( popt[0] )
                list_up.append( popt[3] / popt[2] )
            
        except RuntimeError:
            print 'Optimal parameters not found'

        # Plot data
        plt.plot( xdata[np.where(dx == 'AD')], ydata[np.where(dx == 'AD')], 'o', color='r', label='AD')
        plt.plot( xdata[np.where(dx == 'MCI') ], ydata[np.where(dx == 'MCI')], 'o', color='y', label='MCI')
        plt.plot( xdata[np.where(dx == 'CN')], ydata[np.where(dx == 'CN')], 'o', color='g', label='CN')
        plt.show()
        
    

data_file = os.path.join( adni.project_folder, 'lists/volumes_multiscan.csv' )

vol_index = 75
print 'Analysing', adni.volume_names[vol_index]

baseline_rid = None
last_dx = None
traj_x = []
traj_y = []
traj_dx = []
with open( data_file, 'rb' ) as csvfile:
    reader = csv.reader( csvfile )  
    headers = reader.next()
    for row in reader:
        # Get rid and diagnosis
        rid = int( row.pop(0) )
            
        # Get scan time
        viscode = row.pop(0)
        if viscode == 'bl':
            scan_time = 0
            baseline_rid = rid
        elif re.match('m[0-9][0-9]', viscode):
            scan_time = int(viscode[1:])
        else:
            print 'ERROR: Invalid viscode:', viscode
             
        # Get diagnosis
        dx = row.pop(0)
                      
        # Get diagnosis
        factor = float( row.pop(0) )
        
        # Normalise volumes
        volumes = [float(vol) for vol in row]
#         volumes_sum = sum(volumes)
#         volumes = [vol/volumes_sum for vol in volumes]
#         volumes = [float(vol) * factor for vol in row]
        
        if viscode == 'bl':
            # print 'Plotting subject', rid
            plot_trajectories( traj_x, traj_y, traj_dx )
            traj_x = []
            traj_y = []
            traj_dx = []
        traj_x.append( scan_time )
        traj_y.append( volumes[vol_index] )
        traj_dx.append( dx )

print np.mean( list_r )
print np.mean( list_t0 )
print np.mean( list_up )

x = np.linspace( 0, 48, 100)
y = sigmoid( x, np.mean( list_t0 ), np.exp( np.mean( list_r ) ), 1, np.mean( list_up ) )
plt.plot( x, y, color='r', label='fit1', linewidth=2 )
plt.show()
