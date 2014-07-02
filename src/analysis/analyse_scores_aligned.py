#! /usr/bin/env python
# print __doc__

import os.path
import csv
import re
import sqlite3
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

def normalise_and_append( traj_x, traj_y, traj_d ):   
    ret_x = []
    ret_y = []
    ret_d = []
    if len( traj_x ) > 5:
        
        # Convert and sort data        
        xdata = np.array(traj_x)
        ydata = np.array(traj_y)
        ddata = np.array(traj_d)
        
        args = np.argsort( xdata )
        xdata = xdata[args]
        ydata = ydata[args]
        ddata = ddata[args]
        
#         dydata = []
#         for i in range(len(xdata)-1):
#             dydata.append( float(ydata[i+1] - ydata[i]) / float(xdata[i+1] - xdata[i]) )
            
        if ddata[-1] == 1.0 and ddata[0] == 0.5:
            x_prev = xdata[0];
            for d, x, y in zip(ddata, xdata, ydata):
                if d == 1.0:
                    t_convert = x_prev + (x-x_prev) / 2
                    break
                else:
                    x_prev = x
            
            for x, y, d in zip(xdata, ydata, ddata):
                ret_x.append( x-t_convert )
                ret_y.append( y )
                ret_d.append( d )
#             plt.plot( xdata[1:], dydata, alpha=0.3 )
    return ret_x, ret_y, ret_d


         
#data_file = os.path.join( adni.project_folder, 'lists/volumes_probabilistic.csv' )
data_file = os.path.join( adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv' )

con = sqlite3.connect( os.path.join( adni.project_folder, 'lists', 'adni.db' ) )
con.row_factory = sqlite3.Row
cur = con.cursor()

def analyse_score( name_db, name_hr, est_lower, est_r  ):
    print 'Analysing', name_hr
    list_x = []
    list_y = []
    list_d = []
    traj_x = []
    traj_y = []
    traj_d = []
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
            
            # Get and normalise MMSE
            #cur.execute( "SELECT mmse FROM adnimerge WHERE rid = " + str(rid) + " AND viscode = '" + viscode + "'" )
            cur.execute( "SELECT " + name_db + " FROM adnimerge JOIN adnimerge_cog USING (rid, viscode) WHERE rid = " + str(rid) + " AND viscode = '" + viscode + "'" )
            scans = cur.fetchall()
            if len(scans) == 1:
                try:
                    score = float( scans[0][name_db] )
                except:
                    score = None
    
            if rid != previous_rid and previous_rid != None:
                # Plot previous subject
                traj_x, traj_y, traj_d = normalise_and_append( traj_x, traj_y, traj_d )
                list_x = list_x + traj_x
                list_y = list_y + traj_y
                list_d = list_d + traj_d
                traj_x = []
                traj_y = []
                traj_d = []
            if score != None:
                traj_x.append( scan_time )
                traj_y.append( score )
                traj_d.append( dx )
    
            # Store previous rid for next row
            previous_rid = rid
    
    plt.title( name_hr + ' mean model (' + str(len(list_x)) + ' values)' )
    #plt.title( 'MMSE change' )
    plt.xlabel( 'Months relative to point of convergence' )
    #plt.xlabel( 'Months after baseline' )
    plt.ylabel( name_hr )
    #plt.ylabel( 'MMSE change per month' )
    plt.legend( [mpl.patches.Rectangle((0,0), 1, 1, fc=(0.0, 0.5, 0.0)),
                 mpl.patches.Rectangle((0,0), 1, 1, fc=(0.8, 0.8, 0.0)), 
                 mpl.patches.Rectangle((0,0), 1, 1, fc=(1.0, 0.0, 0.0))],
                ['CN','MCI','AD'],
                bbox_to_anchor=(0.25,0.95) )
    plt.scatter( list_x, list_y, c=list_d, vmin=0.0, vmax=1.0, linewidths=0, cmap=adni.adni_cmap, alpha=0.3 )
    
    # -----------------------
    # Estimate mean model
    # -----------------------
    try:
        # Estimate parameters
        list_x = np.array( list_x )
        list_y = np.array( list_y )
        
        # Fit curve
        popt, _ = curve_fit( exponential, list_x, list_y, p0=(0, est_r, 0.01, est_lower) )
        #popt, _ = curve_fit( sigmoid, list_x, list_y, p0=(0, 1, np.min(list_y), np.max(list_y)) )
        print 'Fitted parameters:', popt
         
        # Plot curve
        x = np.linspace( -50, 50, 200 )
        y = exponential( x, *popt)
        #y = sigmoid( x, *popt)
        plt.plot( x, y )
        #plt.axis( [min(list_x), max(list_x),  min(list_y), max(list_y)] )
    except RuntimeError:
        print 'Optimal parameters not found'             
    
    plt.show()

analyse_score( 'moca', 'MOCA', 30, -20 )
analyse_score( 'mmse', 'MMSE', 30, -20 )
analyse_score( 'cdrsb', 'CDR-SB', 0, 20 )
analyse_score( 'adas11', 'ADAS 11', 0, 20 )
analyse_score( 'adas13', 'ADAS 13', 0, 30 )
analyse_score( 'faq', 'FAQ', 0, 30 )
