import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

(x2,y2,z2) = np.loadtxt("s.dat", unpack=True)

x1 = np.zeros(np.size(x2))
y1 = np.zeros(np.size(y2))
z1 = np.zeros(np.size(z2))

s = np.transpose(np.array([x2,y2,z2]))

#for i in range(0,200):
#	for j in range(0,200):
#		if(i != j and np.sum((s[i]-s[j])**2) < 1.0e-6):x
#			print(i,j, s[i], s[j])
			

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(x2,y2,z2)
#ax.set_xlim([-1, 1])
#ax.set_ylim([-1, 1])
#ax.set_zlim([-1, 1])
plt.show()
