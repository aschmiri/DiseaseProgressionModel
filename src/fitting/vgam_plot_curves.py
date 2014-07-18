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
        
    
def plot( score_name, points_file, curves_file, plot_points, save_file=False, dens_file=None ):
    plt.figure( figsize=(12, 5), dpi=100 )
    
    #
    # Get data
    #
    r = mlab.csv2rec( curves_file )
    r_sorted = np.sort( r, order=r.dtype.names[2] )
    progrs   = r_sorted[r.dtype.names[2]]
    curves   = [r_sorted['x1'],  r_sorted['x5'],  r_sorted['x25'], r_sorted['x50'],
                r_sorted['x75'], r_sorted['x95'], r_sorted['x99']]
    min_progress = np.min(  progrs )
    max_progress = np.max(  progrs )
    progress_offset = (max_progress - min_progress) * 0.1
    
    #
    # Plot densities
    #
    if dens_file != None:
        y, densities = read_dens_file( dens_file )
        
        ax1 = plt.subplot( 1, 2, 1 )
        ax2 = plt.subplot( 1, 2, 2 )
        ax2.set_title( 'Probability density distributions for ' + score_name )
        ax2.set_xlabel( 'Metric value' )
        ax2.set_ylabel( 'Probability' )

            
        min_val = np.min( curves )
        max_val = np.max( curves )
        progr_samples = [-36, -18, 0, 18, 36 ]
        
        sample_cmap = cmx.ScalarMappable(
            norm = colors.Normalize( vmin=-len(progr_samples)+1, vmax=len(progr_samples)-1 ), 
            cmap = plt.get_cmap( adni.adni_cmap ))
        
        for progr in progr_samples:
            sample_color = sample_cmap.to_rgba( progr_samples.index(progr) )
            ax1.axvline( progr, color=sample_color, linestyle='--', alpha=0.8 )
            ax2.set_xlim( min_val, max_val )
            ax2.plot( y, densities[progr + min_progress], label=str(progr), color=sample_color )
            
        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend( handles2, labels2, fontsize=10 )
    else:
        ax1 = plt.subplot( 1, 1, 1 )
    
    #
    # Plot points
    #
    if plot_points:
        m = mlab.csv2rec( points_file )
        progr_points = m['progress']
        value_points = m['value']
        diagn_points = [0.5 if p < 0 else 1.0 for p in progr_points]
        ax1.scatter( progr_points, value_points, c=diagn_points, vmin=0.0, vmax=1.0, linewidths=0, cmap=adni.adni_cmap, alpha=0.25 )
    
    #
    # Plot the percentile curves
    #
    labels   = ['1%', '5%', '25%', '50%', '75%', '95%', '99%']
    greyvals = ['0.45', '0.3', '0.15', '0', '0.15', '0.3', '0.45']
    #styles   = ['g-', 'g-', 'g-', 'b-', 'g-', 'g-', 'g-']

    for (curve, greyval, label) in zip( curves, greyvals, labels ):
        ax1.plot( progrs, curve, color=greyval )
        ax1.text( progrs[-1]+5, curve[-1]-0.5, label, fontsize=10 )
    
    ax1.set_title( 'Percentile curves for ' + score_name )
    ax1.set_xlabel( 'Disease progression relative to point of conversion' )
    ax1.set_ylabel( 'Metric value' )
    ax1.set_xlim( min_progress-progress_offset, max_progress+2*progress_offset )
    
    plt.tight_layout()

    #
    # Draw or save the plot
    #
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
        

