#! /usr/bin/env python
# print __doc__

import os.path
import csv
import numpy as np
import re
from scipy.optimize import curve_fit
import common.adni_tools as adni
import matplotlib.pyplot as plt
import matplotlib as mpl

cdict = {
  'red'  :  ((0.0, 0.0, 0.0), (0.5, 0.8, 0.8), (1.0, 1.0, 1.0)),
  'green':  ((0.0, 0.5, 0.5), (0.5, 0.8, 0.8), (1.0, 0.0, 0.0)),
  'blue' :  ((0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0))
}
my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict)

def sigmoid( t, t0, r, lower, upper ):
    y = lower + (upper - lower) * r * ( 1 / (1 + np.exp(-(t-t0))) )
    return y

def exponential( t, t0, r, lower ):
    y = lower + r * np.exp(t-t0) 
    return y

def plot_trajectories( traj_x, traj_y, traj_d, rid=0, plot_steps=100 ):
    if len( traj_x ) > 5:# and traj_d[-1] == 'CN':
        
        # Convert and sort data        
        xdata = np.array(traj_x)
        ydata = np.array(traj_y)
        ddata = np.array(traj_d)
        
        args = np.argsort( xdata )
        xdata = xdata[args]
        ydata = ydata[args]
        ddata = ddata[args]

        if ddata[-1] == 1.0 and ddata[0] == 0.5:
            x_prev = xdata[0];
            y_prev = ydata[0];
            for d, x, y in zip(ddata, xdata, ydata):
                if d == 1.0:
                    t_convert = x_prev + (x-x_prev) / 2
                    y_convert = y_prev + (y-y_prev) / 2
                    break
                else:
                    x_prev = x
        
            for x, y, d in zip(xdata, ydata, ddata):
                list_x.append( x-t_convert )
                list_y.append( y/y_convert )
                list_d.append( d )
        
#data_file = os.path.join( adni.project_folder, 'lists/volumes_probabilistic.csv' )
data_file = os.path.join( adni.project_folder, 'lists/volumes_segbased.csv' )

for vol_index in range(len(adni.volume_names)): # 18
    print 'Analysing', adni.volume_names[vol_index]
    list_x = []
    list_y = []
    list_d = []
        
    traj_x = []
    traj_y = []
    traj_d = []
    baseline_rid = None
    previous_rid = None
    with open( data_file, 'rb' ) as csvfile:
        rows = csv.DictReader( csvfile )
        for row in rows:
            # Get rid and diagnosis
            rid = int( row['RID'] )
                
            # Get scan time
            viscode = row['VISCODE']
            if viscode == 'bl':
                scan_time = 0
                baseline_rid = rid
            elif re.match('m[0-9][0-9]', viscode):
                scan_time = int( viscode[1:] )
            else:
                print 'ERROR: Invalid viscode:', viscode
                break
                 
            # Get diagnosis
            dx_str = row['DX.scan']
            if dx_str == 'AD':
                dx = 1.0
            elif dx_str == 'MCI':
                dx = 0.5
            elif dx_str == 'CN':
                dx = 0.0
            else:
                print 'ERROR: Invalid diagnosis:', viscode
                break
            
            # Get factor
            factor = 1 / float( row['FactorMNI'] )
            
            # Get and normalise volumes
            volume = float( row[adni.volume_names[vol_index]] )
            # volumes = [float(vol) * factor for vol in row]
    
            if rid != previous_rid and previous_rid != None:
                # Plot previous subject
                print 'Plotting subject', previous_rid, factor
                plot_trajectories( traj_x, traj_y, traj_d, previous_rid )
                traj_x = []
                traj_y = []
                traj_d = []
            traj_x.append( scan_time )
            traj_y.append( volume * factor )
            traj_d.append( dx )
    
            # Store previous rid for next row
            previous_rid = rid
    
    plt.title( adni.volume_names[vol_index] + ' mean model ' )
    plt.xlabel( 'Months after baseline' )
    plt.ylabel( 'Volume growth' )
    plt.legend( [mpl.patches.Rectangle((0,0), 1, 1, fc=(0.0, 0.5, 0.0)),
                 mpl.patches.Rectangle((0,0), 1, 1, fc=(0.8, 0.8, 0.0)), 
                 mpl.patches.Rectangle((0,0), 1, 1, fc=(1.0, 0.0, 0.0))],
                ['CN','MCI','AD'],
                bbox_to_anchor=(0.25,0.95) )
    plt.scatter( list_x, list_y, c=list_d, vmin=0.0, vmax=1.0, linewidths=0, cmap=my_cmap )
    
    
    try:
        # Fit curve
        popt, _ = curve_fit( exponential, list_x, list_y, p0=(0, 0.001, np.mean(list_y)) )
        #popt, _ = curve_fit( sigmoid, list_x, list_y, p0=(0, 1, np.min(list_y), np.max(list_y)) )
        
        # Plot curve
        x = np.linspace( -50, 50, 200 )
        y = exponential( x, *popt)
        #y = sigmoid( x, *popt)
        plt.plot( x, y )
    except RuntimeError:
        print 'Optimal parameters not found'             
    
    
    plt.show()
