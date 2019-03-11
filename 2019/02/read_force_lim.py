import numpy as np
import pandas as pd

lim_path = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z1-001e.lim'

# Load in file.
with open(lim_path) as f:
    lim = f.read()

# This splits at each 'Static forces' instance, where there is one for each
# radial (i.e. ix) location (skip first entry since it's a colon).
lim_sf = lim.split('Static forces')[1:]

# Fix the last element so it doesn't catch everything after ([:-2] to ignore an
# extra \n and space that busg the rest up).
lim_sf[-1] = lim_sf[-1].split('***')[0][:-2]

col_names = ['IX', 'IY', 'XOUT', 'YOUT', 'FEG', 'FIG', 'FF', 'FE',
             'FVH', 'FF2', 'FE2', 'FVH2', 'FTOT1', 'FTOT2', 'TEGS',
             'TIGS', 'CFSS', 'CFVHXS', 'VP1', 'VP2', 'FFB', 'FEB', 'CVHYS',
             'CEYS', 'TE', 'TI', 'NE', 'VELB']

mylist = []
for sf in lim_sf:
    foo = sf.split('\n')[1:]
    foo = [bar.split() for bar in foo]
    df = pd.DataFrame(foo, columns=col_names)
    mylist.append(df)

big_df = pd.concat(mylist, keys=np.arange(1, len(mylist)))
