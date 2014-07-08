#! /usr/bin/env python
# print __doc__

import os.path
import argparse
import csv
import re
import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt
import common.adni_tools as adni

parser = argparse.ArgumentParser()
parser.add_argument( '-s', '--sigmoid', action='store_true', default=False, help='use sigmoid function for fitting' )
a = parser.parse_args()

def sigmoid( t, t0, r, lower, upper ):
    y = lower + (upper - lower) * r * ( 1 / (1 + np.exp(-(t-t0))) )
    return y

def exponential( t, r, r2, lower ):
    y = lower + r * np.exp( r2 * t ) 
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
            value = adni.safe_cast( row[name_csv] )
            if normalise:
                value = value / float( row['FactorMNI'] )
            
            if rid != previous_rid and previous_rid != None:
                # Handle data of previous subject
                if traj_d[-1] == 1.0 and traj_d[0] == 0.5:
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

def compute_binned_data( data_x, data_y, bin_size=6, min_bin_size=10, 
                            return_sample_number=True, remove_outliers=True ):
    bins_x = []
    bins_y = []
    bins_s = []
    bins_n = []

    for curr_x in range( -72, 72, bin_size ):
        bin_x = []
        bin_y = []
        for x, y in zip(data_x, data_y):
            if x >= curr_x and x < curr_x + bin_size:
                bin_x.append( x )
                bin_y.append( y )
        
        if len( bin_x ) > min_bin_size:
            selected_x = np.array( bin_x )
            selected_y = np.array( bin_y )
            
            if remove_outliers:
                mean = np.mean( selected_y )
                sigma = np.std( selected_y )
 
                selected = np.where( np.abs( selected_y - mean ) < 2 * sigma )
                selected_x = selected_x[selected]
                selected_y = selected_y[selected]
                
            bins_x.append( np.mean( selected_x ) )
            bins_y.append( np.mean( selected_y ) )
            bins_s.append( np.std( selected_y ) )
            if return_sample_number:
                bins_n.append( len(selected_x) )

    if return_sample_number:
        return bins_x, bins_y, bins_s, bins_n
    else:
        return bins_x, bins_y, bins_s

def compute_bin_error( data_x, data_y, popt, model=exponential, 
                         bin_size=6, min_bin_size=10 ):
    bins_e = []

    for curr_x in range( -72, 72, bin_size ):
        bin_e = []
        for x, y in zip(data_x, data_y):
            if x >= curr_x and x < curr_x + bin_size:
                y_model = model( x, *popt )
                bin_e.append( np.square( y - y_model ) )
        
        if len( bin_e ) > min_bin_size:
            bins_e.append( np.sqrt( np.mean( np.array( bin_e ) ) ) )
    else:
        return bins_e

def fit_data( data_x, data_y, data_sigma=None, model=exponential ):
    try:
        # Estimate parameters
        data_x = np.array( data_x )
        data_y = np.array( data_y )
            
        # Fit curve
        if data_sigma == None:
            popt, _ = curve_fit( model, data_x, data_y )
        else:
            popt, _ = curve_fit( model, data_x, data_y, sigma=data_sigma )
        print 'Fitted parameters:', popt
    
        return popt

    except RuntimeError:
        print 'Optimal parameters not found'
        return None

def evaluate_robustness( data_x, data_sigma, popt, model=exponential ):
    delta = 0.0001
    data_robustness = []
    for x, s in zip(data_x, data_sigma):
        # Approximate slope with finite differences
        y2 = (model( x+delta, *popt ) - model( x-delta, *popt )) / (2*delta)
        data_robustness.append( np.abs( y2 ) / s )
    
    return np.array( data_robustness )

def evaluate_group_separation( data_x, data_y, popt, model=exponential ):
    data_x = np.array( data_x )
    data_y = np.array( data_y )
    
    y_mci = data_y[ np.where( data_x < 0 ) ]
    y_ad  = data_y[ np.where( data_x > 0 ) ]
    y_conv = model( 0, *popt )
    mean_mci = np.mean( y_mci )
    mean_ad  = np.mean( y_ad )
    
    if mean_mci < mean_ad:
        err_mci = len( np.where( y_mci > y_conv )[0] )
        err_ad  = len( np.where( y_ad  < y_conv )[0] )
    else:
        err_mci = len( np.where( y_mci < y_conv )[0] )
        err_ad  = len( np.where( y_ad  > y_conv )[0] )
    return float(err_mci + err_ad) / float(len(data_y)) 
        
def analyse_metric( name_csv, name_hr, normalise, plot=False ):
    print 'Analysing', name_hr
    collect_data( name_csv, normalise )
    age_regress_data()

    global list_x
    global list_y
    global list_d
    global metric_name
    global metric_robustness_min
    global metric_robustness_max
    global metric_robustness_mean
    global metric_group_seperation

    if plot:
        _, ax1 = plt.subplots()
        plt.title( name_hr + ' mean model (' + str(len(list_x)) + ' values)' )
        ax1.set_xlabel( 'Months relative to point of conversion' )
        ax1.set_ylabel( name_hr )
        ax1.legend( [mpl.patches.Rectangle((0,0), 1, 1, fc=(0.0, 0.5, 0.0)),
                     mpl.patches.Rectangle((0,0), 1, 1, fc=(0.8, 0.8, 0.0)), 
                     mpl.patches.Rectangle((0,0), 1, 1, fc=(1.0, 0.0, 0.0))],
                    ['CN','MCI','AD'],
                    bbox_to_anchor=(0.25,0.95) )
        ax1.scatter( list_x, list_y, c=list_d, vmin=0.0, vmax=1.0, linewidths=0, cmap=adni.adni_cmap, alpha=0.25 )
    
    # -----------------------
    # Compute and plot binned fit
    # -----------------------
    bins_x, bins_y, bins_s, bins_n = compute_binned_data( list_x, list_y, bin_size=3,
                                                          return_sample_number=True, 
                                                          remove_outliers=True )
    if plot:
        ax1.errorbar( bins_x, bins_y, yerr=bins_s, 
                      linewidth=0, elinewidth=1, marker='x', c='k', ecolor='0.5' )

    # Scale bins with number of samples
    bins_weight =  [ bins_s[i] / np.sqrt( bins_n[i] ) for i in range(len(bins_s)) ]
    
    # Fit bin data
    popt = fit_data( bins_x, bins_y, data_sigma=bins_weight, model=model )
    if popt != None:
        if plot:
            # Plot fitted curve
            x = np.linspace( -50, 50, 200 )
            y = model( x, *popt)
            ax1.plot( x, y, color='k' )
        
        # Estimate and plot robustness
        bins_e = compute_bin_error( list_x, list_y, popt, model=model, bin_size=3 )
        bins_r = evaluate_robustness( bins_x, bins_e, popt, model=model )
        group_sep = evaluate_group_separation( list_x, list_y, popt, model=model )

        print 'Robustness min:  ', np.min( bins_r )
        print 'Robustness max:  ', np.max( bins_r )
        print 'Robustness mean: ', np.mean( bins_r )
        print 'Group separation:', group_sep
                
        if plot:
            # Plot mean squared error
#             bins_y_model = [model( x, *popt ) for x in bins_x ]
#             ax1.errorbar( bins_x, bins_y_model, yerr=bins_e, 
#                           linewidth=0, elinewidth=1, marker='x', c='r', ecolor='r' )
            
            # Plot robustness
            ax2 = plt.twinx()
            ax2.set_ylim( ymin=0, ymax=0.2 )
            ax2.set_ylabel( 'Robustness' )
            ax2.plot( bins_x, bins_r, color='g' )
        else:
            metric_name.append( name_csv )
            metric_robustness_min.append(  np.min( bins_r ) )
            metric_robustness_max.append(  np.max( bins_r ) )
            metric_robustness_mean.append( np.mean( bins_r ) )
            metric_group_seperation.append( group_sep )

    if plot:
#         ymin = 1 + 1.1 * ( min(list_y) - 1 )
#         ymax = 1 + 1.1 * ( max(list_y) - 1 )
#         ax1.set_ylim( ymin=ymin, ymax=ymax )
    
        plt.show()


#
# Determine fitting model
#
if a.sigmoid:
    model = sigmoid
else:
    model = exponential

#
# 
#
metric_name             = []
metric_robustness_min   = []
metric_robustness_max   = []
metric_robustness_mean  = []
metric_group_seperation = []

analyse_metric( 'MMSE', 'MMSE', False )
analyse_metric( 'CDRSB', 'CDR-SB', False )
analyse_metric( 'ADAS11', 'ADAS 11', False )
analyse_metric( 'ADAS13', 'ADAS 13', False )
analyse_metric( 'FAQ', 'FAQ', False )
for vol_index in range(len(adni.volume_names)):
    analyse_metric( adni.volume_names[vol_index], adni.volume_names[vol_index], True )

metric_name             = np.array( metric_name )
metric_robustness_min   = np.array( metric_robustness_min )
metric_robustness_max   = np.array( metric_robustness_max )
metric_robustness_mean  = np.array( metric_robustness_mean )
metric_group_seperation = np.array( metric_group_seperation )

args = np.argsort( metric_robustness_min )[::-1]
sorted_metric_names = metric_name[args]
print sorted_metric_names
for name in sorted_metric_names:
    if name in ['MMSE', 'CDRSB', 'ADAS11', 'ADAS13', 'FAQ' ]:
        normalise = False
    else:
        normalise = True
    analyse_metric( name, name, normalise, plot=True )


