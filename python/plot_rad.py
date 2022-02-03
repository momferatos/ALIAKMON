import os
import sys
import numpy as np
import matplotlib.pyplot as plt

dirs = sys.argv[1:]
labels = dirs

fig, axs = plt.subplots(1, 2)
axs[0].set_title('Radiative heat flux integrated on a sphere')
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('$Q_r$ (W)')
axs[1].set_title('Total radiative energy')
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('$E_r$ (J)')
for idir, pdir in enumerate(dirs):
    x, y, z = np.loadtxt(os.path.join(pdir, 'fort.432'), unpack=True)
    axs[0].plot(x, y, label=labels[idir])
    axs[1].plot(x, z, label=labels[idir])
plt.legend(loc=0)
plt.show()
