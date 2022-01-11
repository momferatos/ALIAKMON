import numpy as np
import matplotlib.pyplot as plt

x, y, z = np.loadtxt('fort.432', unpack=True)
fig, axs = plt.subplots(1, 2)
axs[0].plot(x, y)
axs[1].plot(x, z)
plt.show()
