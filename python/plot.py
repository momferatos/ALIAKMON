import numpy as np
from vtk import vtkStructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys
import glob

minval = 1.0e-3

ncmap = np.load('/users/pr008/gmomfer/codes/aliakmon/python/turb_cmap.npy') / 255.
cmap = ListedColormap(ncmap)

reader = vtkStructuredGridReader()

if len(sys.argv) == 3:
    vtkfile = sys.argv[2]
else:
    vtkfiles = glob.glob('./*.vtk', recursive=False)
    vtkfiles.sort(reverse=True)
    vtkfile=vtkfiles[0]

reader.SetFileName(vtkfile)
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()

fieldkey = sys.argv[1]

    

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

field = VN.vtk_to_numpy(data.GetPointData().GetArray(fieldkey))    
field = field.reshape(dim)[:,:,0]
field = np.where(field < 0.0001, 0.0001, field)



#rcParams.update({'figure.figsize':figsize})
fig = plt.figure()
plt.axes([0,0,1,1]) # Make the plot occupy the whole canvas
plt.axis('off')
#fig.set_size_inches(figsize)
plt.contourf(np.log(field), 256, cmap=cmap)
plt.axis('equal')
plt.show()
    
