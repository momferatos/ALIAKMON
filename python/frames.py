#!/usr/bin/env python3
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys
import glob
import os

minval = 1.0e-3

pydir = os.path.dirname(os.path.realpath(__file__))

ncmap_blueblack = np.load(pydir + '/blue-black_cmap.npy') / 255.
ncmap_blackbody = np.load(pydir + '/black-body_cmap.npy') / 255.

cmap_blueblack = ListedColormap(ncmap_blueblack)
cmap_blackbody = ListedColormap(ncmap_blackbody)

paths = sys.argv[1:]

nh5files = 0
for ipath,path in enumerate(paths):
    h5files = glob.glob(path + '/*.h5', recursive=False)
    nh5files = nh5files + len(h5files)


ih5file = 0
for ipath,path in enumerate(paths):
    h5files = glob.glob(path + '/*.h5', recursive=False)
    nh5files = len(h5files)
    
    h5files.sort()
    maxfield = {}
    minfield = {}
    for ih5file,h5file in enumerate(h5files):
        with h5py.File(h5file, 'r') as h5file:
            (print('Calculating maxima/minima: # '
                   + f'{ih5file+1:d}/{nh5files:d}', end='\r'))
            ih5file += 1
            fieldkeys = h5file.keys() 
            fieldkeys = ([fieldkey for fieldkey in fieldkeys
                          if fieldkey != 'time'])
            for fieldkey in fieldkeys:
                field = np.array(h5file[fieldkey])
                min_value = np.finfo(field.dtype).min
                max_value = np.finfo(field.dtype).max
                maxfield[fieldkey] = max(maxfield.setdefault(fieldkey,
                                        min_value), np.max(field))
                minfield[fieldkey] = min(minfield.setdefault(fieldkey,
                                        max_value), np.min(field))
    print('Done calculating maxima/minima.                                   ')
    pngdir = os.path.join(path, 'png')
    if not os.path.isdir(pngdir):
        os.mkdir(pngdir)
    ih5file = 0
    for ih5file,h5file in enumerate(h5files):        
        with h5py.File(h5file, 'r') as h5file:            
            print('Processing file # ' + f'{ih5file+1:d}/{nh5files:d}',
                  end='\r')
            for fieldkey in fieldkeys:
                pngfile = os.path.join(pngdir, f'{fieldkey}.{ih5file:06}.png')
                if os.path.isfile(pngfile):
                    continue
                field = np.array(h5file[fieldkey])
                field2 = field[:,0:1920-1024]
                field = np.concatenate((field,field2), axis=1)
                field2 = field[0:1080-1024,:]
                field = np.concatenate((field,field2), axis=0)
                scale = (maxfield[fieldkey] - minfield[fieldkey]) ** (-1)
                field = scale * (field - minfield[fieldkey])
                if 'scl' in fieldkey:
                    cmap = cmap_blackbody
                else:
                    cmap = cmap_blueblack
                    field = np.log(np.where(field < minval, minval, field))
                #cmap = cmap_blueblack
#                cmap = cmap_blackbody if 'scl' in fieldkey else cmap_blueblack
                plt.imsave(pngfile, field, cmap=cmap)
    print(f'Finished directory {path}.                                      ')


