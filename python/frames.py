import numpy as np
from vtk import vtkStructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys
import glob

minval = 1.0e-3

ncmap = np.load('./turb_cmap.npy') / 255.

cmap = ListedColormap(ncmap)

reader = vtkStructuredGridReader()

paths = sys.argv[1:]

nvtkfiles = 0
for ipath,path in enumerate(paths):
    vtkfiles = glob.glob(path + '/*.vtk', recursive=False)
    nvtkfiles = nvtkfiles + len(vtkfiles)


ivtkfile = 0
for ipath,path in enumerate(paths):
    vtkfiles = glob.glob(path + '/*.vtk', recursive=False)
    nvtkfiles = nvtkfiles + len(vtkfiles)
    
    vtkfiles.sort()
  
    for vtkfile in vtkfiles:
        
        reader.SetFileName(vtkfile)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()

        data = reader.GetOutput()
        
        print('Calculating maxima: # ' + f'{ivtkfile+1:d}/{nvtkfiles:d}', end='\r')
        ivtkfile += 1
        fieldkeys = ['w', 'e']
        maxfield = {}
        for fieldkey in fieldkeys:
            field = VN.vtk_to_numpy(data.GetPointData().GetArray(fieldkey))
            maxfield[fieldkey] = np.max(field)

    ivtkfile = 0
    for vtkfile in vtkfiles:
        tmp = vtkfile.split('-')
        tmp = tmp[-1].split('.')

        reader.SetFileName(vtkfile)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()

        data = reader.GetOutput()

        dim = data.GetDimensions()
        vec = list(dim)
        scl = vec
        vec = [i-1 for i in dim]
        vec.append(3)
        print('Processing file # ' + f'{ivtkfile+1:d}/{nvtkfiles:d}', end='\r')
        for fieldkey in fieldkeys:
            field = VN.vtk_to_numpy(data.GetPointData().GetArray(fieldkey))    
            field = field.reshape(dim)[:,:,0]
            # field2 = field[:,0:1920-1024]
            # field = np.concatenate((field,field2), axis=1)
            # field2 = field[0:1080-1024,:]
            # field = np.concatenate((field,field2), axis=0)
            field = field / maxfield[fieldkey]
            field = np.log(np.where(field < minval, minval, field))
            plt.imsave(fieldkey + f'.{ivtkfile:06}.png', field, cmap=cmap)
        ivtkfile += 1


ivtkfile=0
for ipath,path in enumerate(paths):
    vtkfiles = glob.glob(path + '/*.vtk', recursive=False)

    vtkfiles.sort()

    for vtkfile in vtkfiles:

        reader.SetFileName(vtkfile)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        print('Calculating maxima: # ' + f'{ivtkfile+1:d}/{nvtkfiles:d}', end='\r')
        fieldkeys = ['w', 'e']
        maxfield = {'w':0.0, 'e':0.0}
        for fieldkey in fieldkeys:
            data = reader.GetOutput()
            field = VN.vtk_to_numpy(data.GetPointData().GetArray(fieldkey))
            maxfield[fieldkey] = max(maxfield[fieldkey], np.max(field))


    for vtkfile in vtkfiles:
        
        reader.SetFileName(vtkfile)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()

        data = reader.GetOutput()

        dim = data.GetDimensions()
        vec = list(dim)
        scl = vec
        vec = [i-1 for i in dim]
        vec.append(3)
        print('Processing files: # ' + f'{ivtkfile+1:d}/{nvtkfiles:d}', end='\r')
        for fieldkey in fieldkeys:
            field = VN.vtk_to_numpy(data.GetPointData().GetArray(fieldkey))    
            field = field.reshape(dim)[:,:,0]
            # field2 = field[:,0:1920-1024]
            # field = np.concatenate((field,field2), axis=1)
            # field2 = field[0:1080-1024,:]
            # field = np.concatenate((field,field2), axis=0)
            field = field / maxfield[fieldkey]
            field = np.log(np.where(field < minval, minval, field))
            plt.imsave(fieldkey + f'.{ivtkfile:06}.png', field, cmap=cmap)
        ivtkfile += 1
