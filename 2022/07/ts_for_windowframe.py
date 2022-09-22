# Script to give some copy/paste output for window_frame.xlsx so we can compare
# the TS data to the RCP.
import pickle
import numpy as np


tspath = "/Users/zamperini/Documents/d3d_work/divimp_files/ts_190442.pickle"
with open(tspath, "rb") as f:
    ts = pickle.load(f)

time_window = [4390, 4550]
time = ts["core"]["time"]
mask = np.logical_and(time>=time_window[0], time<=time_window[1])
te = ts["core"]["te"][:, mask].flatten()
ne = ts["core"]["ne"][:, mask].flatten()
psin = ts["core"]["psin"][:, mask].flatten()
