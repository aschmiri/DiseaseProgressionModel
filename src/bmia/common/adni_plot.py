#! /usr/bin/env python2.7
import matplotlib as mpl
progression_dict = {
    'red': ((0.0, 0.0, 0.0), (0.5, 0.8, 0.8), (1.0, 1.0, 1.0)),
    'green': ((0.0, 0.5, 0.5), (0.5, 0.8, 0.8), (1.0, 0.0, 0.0)),
    'blue': ((0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0))
}
progression_cmap = mpl.colors.LinearSegmentedColormap('my_colormap', progression_dict)
