#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys
import os
import h5py

minval = 1.0e-3

pydir = os.path.dirname(os.path.realpath(__file__))

ncmap_blueblack = np.load(pydir + '/blue-black_cmap.npy') / 255.
ncmap_blackbody = np.load(pydir + '/black-body_cmap.npy') / 255.

cmap_blueblack = ListedColormap(ncmap_blueblack)
cmap_blackbody = ListedColormap(ncmap_blackbody)

fname = sys.argv[1]
h5key = sys.argv[2]

if not os.path.isfile(fname):
    print(f'error: can\'t open {fname}.')
    sys.exit(1)
    
with h5py.File(fname, 'r') as h5file:
    h5keys = h5file.keys()
    if h5key not in h5keys:
        print(f'error: key {h5key} not found')
        sys.exit(1)
    field = np.array(h5file[h5key])
    if 'scl' in h5key or h5key == 'G':
        cmap = cmap_blackbody
    else:
        cmap = cmap_blueblack
        field = np.log(np.where(field < minval, minval, field))
    plt.figure()
    plt.axes([0,0,1,1]) # Make the plot occupy the whole canvas
    plt.axis('off')
    plt.contourf(field, 256, cmap=cmap)
    plt.axis('equal')
    plt.show()
    
