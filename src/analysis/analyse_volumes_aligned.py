#! /usr/bin/env python
# print __doc__

import os.path
import csv
import re
import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt
import common.adni_tools as adni

def sigmoid( t, t0, r, lower, upper ):
    y = lower + (upper - lower) * r * ( 1 / (1 + np.exp(-(t-t0))) )
    return y

def exponential( t, t0, r, r2, lower ):
    y = lower + r * np.exp( r2 * t - t0 ) 
    return y

def normalise_and_append( traj_x, traj_y, traj_d, traj_a, normalise ):
    global list_x
    global list_y
    global list_d
    global list_a

    # Convert and sort data        
    xdata = np.array(traj_x)
    ydata = np.array(traj_y)
    ddata = np.array(traj_d)
    adata = np.array(traj_a)
    
    args = np.argsort( xdata )
    xdata = xdata[args]
    ydata = ydata[args]
    ddata = ddata[args]
    adata = adata[args]

    x_prev = xdata[0]
    y_prev = ydata[0]
    for d, x, y in zip(ddata, xdata, ydata):
        if d == 1.0:
            t_convert = x_prev + (x-x_prev) / 2
            y_convert = y_prev + (y-y_prev) / 2
            break
        else:
            x_prev = x

    if not normalise:
        y_convert = 1

    for x, y, d, a in zip(xdata, ydata, ddata, adata):
        list_x.append( x-t_convert )
        list_y.append( y/y_convert )
        list_d.append( d )
        list_a.append( a )

def collect_data( name_csv, normalise ):
    global list_x
    global list_y
    global list_d
    global list_a

    list_x = []
    list_y = []
    list_d = []
    list_a = []
    
    traj_x = []
    traj_y = []
    traj_d = []
    traj_a = []
    
    previous_rid = None
    data_file = os.path.join( adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv' )
    with open( data_file, 'rb' ) as csvfile:
        rows = csv.DictReader( csvfile )
        for row in rows:
            # Get rid and age
            rid = int( row['RID'] )
            age = float( row['AGE.scan'] )
            
            # Get scan time
            viscode = row['VISCODE']
            if viscode == 'bl':
                scan_time = 0
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
            
            # Get and normalise volumes
            value = adni.safe_cast_to_float( row[name_csv] )
            if normalise:
                value = value / float( row['FactorMNI'] )
            
            if rid != previous_rid and previous_rid != None:
                # Plot previous subject
                if traj_d[-1] == 1.0 and traj_d[0] == 0.5:
                    print 'Plotting subject', previous_rid
                    normalise_and_append( traj_x, traj_y, traj_d, traj_a, normalise )
                traj_x = []
                traj_y = []
                traj_d = []
                traj_a = []
            if value != None:
                traj_x.append( scan_time )
                traj_y.append( value )
                traj_d.append( dx )
                traj_a.append( age )
    
            # Store previous rid for next row
            previous_rid = rid
            
def age_regress_data():
    global list_y
    global list_a
    
    import scipy.stats as stats
    list_y = np.array( list_y )
    list_a = np.array( list_a )
    mean_y = np.mean( list_y )
    
    slope, intercept, _, _, _ = stats.linregress( list_a , list_y )
    list_y  = list_y - (list_a * slope + intercept) + mean_y

def compute_binned_data( data_x, data_y, bin_size=6, min_bin_size=10, scale_with_samples=True ):
    bins_x = []
    bins_y = []
    bins_s = []

    for curr_x in range( -72, 72, bin_size ):
        bin_x = []
        bin_y = []
        for x, y in zip(data_x, data_y):
 
            if x >= curr_x and x < curr_x + bin_size:
                bin_x.append( x )
                bin_y.append( y )
        if len( bin_x ) > min_bin_size:    
            bins_x.append( np.mean( bin_x ) )
            bins_y.append( np.mean( bin_y ) )
            s = np.std( bin_y )
            if scale_with_samples:
                s = s / np.sqrt(len(bin_x))
            bins_s.append( s )
    return bins_x, bins_y, bins_s
    
def fit_data( data_x, data_y, data_sigma=None, model=exponential ):
    try:
        # Estimate parameters
        data_x = np.array( data_x )
        data_y = np.array( data_y )
    
        mci_mean = np.mean( data_y[np.where(data_x<0)] )
        ad_mean  = np.mean( data_y[np.where(data_x>0)] )
        
        print mci_mean, ad_mean
        
        # Fit curve
        if data_sigma == None:
            popt, _ = curve_fit( model, data_x, data_y )
        else:
            popt, _ = curve_fit( model, data_x, data_y, sigma=data_sigma )
        
        print 'Fitted parameters:', popt

        # Plot curve
        x = np.linspace( -50, 50, 200 )
        y = model( x, *popt)

        return x, y

    except RuntimeError:
        print 'Optimal parameters not found'
        return None, None

def analyse_metric( name_csv, name_hr, normalise ):
    print 'Analysing', name_hr
    collect_data( name_csv, normalise )
    #age_regress_data()

    global list_x
    global list_y
    global list_d
    
    plt.title( name_hr + ' mean model (' + str(len(list_x)) + ' values)' )
    plt.xlabel( 'Months relative to point of conversion' )
    plt.ylabel( name_hr )
    plt.legend( [mpl.patches.Rectangle((0,0), 1, 1, fc=(0.0, 0.5, 0.0)),
                 mpl.patches.Rectangle((0,0), 1, 1, fc=(0.8, 0.8, 0.0)), 
                 mpl.patches.Rectangle((0,0), 1, 1, fc=(1.0, 0.0, 0.0))],
                ['CN','MCI','AD'],
                bbox_to_anchor=(0.25,0.95) )
    plt.scatter( list_x, list_y, c=list_d, vmin=0.0, vmax=1.0, linewidths=0, cmap=adni.adni_cmap, alpha=0.3 )
    
    # -----------------------
    # Compute bin data
    # -----------------------
    bins_x, bins_y, bins_s = compute_binned_data( list_x, list_y, bin_size=3, scale_with_samples=True )
    plt.errorbar(bins_x, bins_y, yerr=bins_s, linewidth=0, elinewidth=1, marker='x', c='k', ecolor='0.5' )
    
    x, y = fit_data( bins_x, bins_y, bins_s, model=exponential )
    if x != None and y != None:
        plt.plot( x, y, color='k' )
    
    x, y = fit_data( list_x, list_y, model=exponential )
    if x != None and y != None:
        plt.plot( x, y, '--', color='0.5' )

    plt.show()

# analyse_metric( 'MMSE', 'MMSE', False )#30, -20 )
# analyse_metric( 'CDRSB', 'CDR-SB', False )# 0, 20 )
# analyse_metric( 'ADAS11', 'ADAS 11', False )# 0, 20 )
# analyse_metric( 'ADAS13', 'ADAS 13', False )# 0, 30 )
# analyse_metric( 'FAQ', 'FAQ', False )# 0, 30 )
for vol_index in range(len(adni.volume_names)): # 18
    analyse_metric( adni.volume_names[vol_index], adni.volume_names[vol_index], True )
