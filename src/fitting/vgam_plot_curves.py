#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cmx
import common.adni_tools as adni


def read_dens_file( dens_file ):
    import csv
    densities = []
    with open( dens_file, 'rb' ) as csvfile:
        rows = csv.reader(csvfile)
        rows.next()
        for row in rows:
            name = row[0]
            data = [ float(row[i]) for i in range(1,len(row))]
            if name == 'values':
                y = data
            else:
                densities.append( data )
    return np.asarray( y ), np.asarray( densities )
        
    
def plot( score_name, points_file, curves_file, plot_points, save_file=False, dens_file=None, min_progression=-42 ):
    # Load csv in a structured array
    r = mlab.csv2rec( curves_file )
    plt.figure( figsize=(13, 7), dpi=100 )
    
    # Sort the array on age
    r_sorted = np.sort( r, order=r.dtype.names[2] )
    progrs   = r_sorted[r.dtype.names[2]]
    curves   = [r_sorted['x1'],  r_sorted['x5'],  r_sorted['x25'], r_sorted['x50'],
                r_sorted['x75'], r_sorted['x95'], r_sorted['x99']]
    labels   = ['1%', '5%', '25%', '50%', '75%', '95%', '99%']
    greyvals = ['0.6', '0.4', '0.2', '0', '0.2', '0.4', '0.6']
    x_offset = 0.5
    y_offset = 0.5
    
    if dens_file != None:
        y, densities = read_dens_file( dens_file )
        
        ax1 = plt.subplot(1,2,1)
        ax2 = plt.subplot(1,2,2)
        ax2.set_title( 'Probability density distributions for ' + score_name )
        ax2.set_xlabel( 'Metric value' )
        ax2.set_ylabel( 'Probability' )
            
        min_val = np.min(curves)
        max_val = np.max(curves)
        progr_samples = [-36, -18, 0, 18, 36 ]
        print min_val,max_val
        
        sample_cmap = cmx.ScalarMappable(
            norm = colors.Normalize( vmin=0, vmax=len(progr_samples)-1 ), 
            cmap = plt.get_cmap('winter') )
        
        for progr in progr_samples:
            sample_color = sample_cmap.to_rgba( progr_samples.index(progr) )
            ax1.axvline( progr, color=sample_color, linestyle='--', alpha=0.5 )
            ax2.set_xlim( min_val, max_val )
            ax2.plot( y, densities[progr + min_progression], label=str(progr), color=sample_color )
    else:
        ax1 = plt.subplot(111)
    
    if plot_points:
        # Collect all measurement points from the measurement file
        m = mlab.csv2rec( points_file )
        progr_points = m['progress']
        value_points = m['value']
        ax1.scatter( progr_points, value_points, vmin=0.0, vmax=1.0, linewidths=0, cmap=adni.adni_cmap, alpha=0.25 )
    
    # Plot the percentile curves
    for (curve, greyval, label) in zip( curves, greyvals, labels ):
        ax1.plot( progrs, curve, color=greyval )
        ax1.text( progrs[-1]+x_offset, curve[-1]-y_offset, label )
    
    ax1.set_title( 'Percentile curves for ' + score_name )
    ax1.set_xlabel( 'Disease progression relative to point of conversion' )
    ax1.set_ylabel( 'Metric value' )
    plt.tight_layout()

    # Draw the plot
    if save_file:
        plot_filename = curves_file.replace( '.csv', '_plot.png' )
        plt.savefig( plot_filename, dpi=100 )
    else:  
        plt.show()

if __name__ == "__main__":
    plot_points = True
    
    for score_name in adni.biomarker_names:
        print 'Generating plot for', score_name
        
        points_file = os.path.join( adni.project_folder, 'data', score_name.replace( ' ', '_' ) + '.csv' )
        curves_file = points_file.replace( '.csv', '_curves.csv' )
        dens_file = points_file.replace( '.csv', '_densities.csv' )
        if os.path.isfile( points_file ) and os.path.isfile( curves_file ): 
            plot( score_name, points_file, curves_file, plot_points, dens_file=dens_file )
        

