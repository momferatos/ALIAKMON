import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

path = sys.argv[1]

series = pd.read_table(path + '/' + 'hydro.dat', delimiter='|', skiprows=1)

series.plot(x='t',y='e')
#series.plot(x='t',y='ke')
plt.show()
# hydro_dat = open(path + '/' + 'hydro.dat', 'r')

# descriptions = hydro_dat.readline().split(separator='|').strip()
# keys = hydro_dat.readline().split(separator='|').strip()

# print(descriptions)
# print(keys)
