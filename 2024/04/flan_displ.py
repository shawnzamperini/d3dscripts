import pickle
import matplotlib.pyplot as plt
import numpy as np


# Load displacement dictionary. 
path = "/Users/zamperini/flandir/reg_testcase1/saved_results_7/displ_dict.pickle"
with open(path, "rb") as f:
    d = pickle.load(f)

# Load lists as numpy arrays.
xstart = np.array(d["xstarts"])
zstart = np.array(d["zstarts"])
xend = np.array(d["xends"])
zend = np.array(d["zends"])

# Displacements.
dx = np.abs(xend - xstart)
dz = np.abs(zend - zstart)

# Probably want to do a spline fit. 

fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.scatter(zstart, dx/dz)
fig.tight_layout()
fig.show()
