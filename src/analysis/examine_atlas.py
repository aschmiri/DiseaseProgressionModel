#! /usr/bin/env python
# print __doc__

import os.path
import argparse
import numpy as np
import nibabel as nib
import matplotlib.pyplot as pyplot
import matplotlib.cm as cm

parser = argparse.ArgumentParser()
parser.add_argument( '-i', '--iteration', dest='iteration', type=int, default=1 )
parser.add_argument( '--min', dest='state_min', type=float, default=0 )
parser.add_argument( '--max', dest='state_max', type=float, default=10 )
parser.add_argument( '--steps', dest = 'state_steps', type=int, default=11 )
a = parser.parse_args()
    
#for mmse in range(65, 85):
class EventHandler:
    def __init__(self, subplot, atlases):
        self.ax = subplot
        self.state = 0
        self.slice = 90
        
        self.update()
            
    def onkeypress(self, event):
        key = event.key
        if key == 'q':
            self.slice = np.clip(self.slice+1, 0, self.slices-1)
        elif key == 'w':
            self.slice = np.clip(self.slice-1, 0, self.slices-1)
        self.update()  
            
    def onscroll(self, event):
        if event.button == 'up':
            self.state = np.clip(self.state-1, 0, a.state_steps-1)
        else:
            self.state = np.clip(self.state+1, 0, a.state_steps-1)
        self.update()

    def update(self):
        self.ax.cla()
        self.ax.set_ylabel('slice %s' % self.slice)
        self.ax.set_xlabel('current state: %s' % self.state)
        
        self.imagedata = nib.load( atlases[self.state] ).get_data()        
        self.im = self.ax.imshow(self.imagedata[:,:,self.slice], interpolation = 'nearest', cmap = cm.gray ) #
        self.im.axes.figure.canvas.draw()
        
atlas_base_folder = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas'
atlases = []
for state in np.linspace( a.state_min, a.state_max, a.state_steps ):
    atlas_folder =  os.path.join( atlas_base_folder, 'model_0' )
    #atlas = os.path.join( atlas_folder, 'average_image_transformed.nii.gz' )
    #state = '{0}'.format(str(round(state,2) if state % 1  else int(state)),1)
    atlas = os.path.join( atlas_folder, 'atlas_m24_AD_' + str(state) + '_average_image_transformed.nii.gz' )
    atlases.append( atlas )

###################################
# Plot image and features
fig = pyplot.figure()
subplot = fig.add_subplot(1,1,1)
tracker = EventHandler( subplot, atlases )

fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
fig.canvas.mpl_connect('key_press_event', tracker.onkeypress)
pyplot.show()


