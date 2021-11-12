#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

fname = sys.argv[1]
if not os.path.isfile(fname):
    print(f'error: can\'t open {fname}.')
    sys.exit(1)
    
ys = sys.argv[2:]

df = pd.read_table(fname, delimiter='|', skiprows=1)

with open(fname, 'r') as hydro_dat:
    descs = list(hydro_dat)[0].split('|')
    
if not ys:
    cols = list(df.columns)
    for col,desc in zip(cols, descs):
        print(f'{col.strip():<10} {desc.strip():>40}')
else:
    for y in ys:
        df.plot(x='t', y=y)
        
#series.plot(x='t',y='ke')
plt.show()
# hydro_dat = open(path + '/' + 'hydro.dat', 'r')

# descriptions = hydro_dat.readline().split(separator='|').strip()
# keys = hydro_dat.readline().split(separator='|').strip()

# print(descriptions)
# print(keys)
