#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import common.adni_tools as adni

def plot( score_name, points_file, curves_file, plot_points ):
    # Load csv in a structured array
    r = mlab.csv2rec( curves_file )
    plot_filename = curves_file.replace( '.csv', '_plot.png' ) 
    plt.figure( figsize=(13, 7), dpi=100 )
    
    if plot_points:
        # Collect all measurement points from the measurement file
        m = mlab.csv2rec( points_file )
        age_points = m['scan_age']
        vol_points = m['volume']
        #plt.plot( age_points, vol_points, 'mo', mfc='none', markersize=4 )
        plt.scatter( age_points, vol_points, vmin=0.0, vmax=1.0, linewidths=0, cmap=adni.adni_cmap, alpha=0.25 )
        
    # Sort the array on age
    r_sorted = np.sort(r, order=r.dtype.names[2])

    age = r_sorted[r.dtype.names[2]]
    x5  = r_sorted['x5']
    x25 = r_sorted['x25']
    x50 = r_sorted['x50']
    x75 = r_sorted['x75']
    x95 = r_sorted['x95']
    
    curves     = [x5,   x25,   x50,   x75,   x95]
    linestyles = ['g-', 'g-',  'b-',  'g-',  'g-']
    labels     = ['5%', '25%', '50%', '75%', '95%']
    x_offset = 0.5
    y_offset = 0.5

    # Plot the percintile curves
    for ( curve, linestyle, label ) in zip( curves, linestyles, labels ):
        plt.plot( age, curve, linestyle, linewidth=2.0)
        plt.text( age[-1] + x_offset, curve[-1] - y_offset, label )
    plt.title( 'Percentile curves for ' + score_name )
    plt.xlabel( 'Time relative to point of conversion' )
    plt.ylabel( 'Volume' )
    plt.grid()
    plt.tight_layout()

    # Draw the plot
#     plt.show()
    plt.savefig( plot_filename, dpi=100 )


if __name__ == "__main__":
    plot_points = True
    
    for score_name in adni.biomarker_names:
        print 'Generating plot for', score_name
        
        points_file = os.path.join( adni.project_folder, 'data', score_name.replace( ' ', '_' ) + '.csv' )
        curves_file = points_file.replace( '.csv', '_curves.csv' )
        if os.path.isfile( points_file ) and os.path.isfile( curves_file ): 
            plot( score_name, points_file, curves_file, plot_points )
        

