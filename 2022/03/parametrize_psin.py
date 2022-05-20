# This script takes in a gfile, and will parameterize a flux tube at a given
# psin into (R, Z) coordinates. This is then to be passed into MAFOT to get
# calculations of the connection lengths.
import pickle
import numpy as np
from LimWallToolkit import LimWallToolkit
import matplotlib.pyplot as plt
from tqdm import tqdm


# Inputs.
#psins = np.linspace(1.13, 1.40, 50)   # 167196
psins = np.linspace(1.18, 1.45, 50)   # 186914
shot = 176971
pol_lims = True

# Choose correct paths for everything. There actually is no reason to separate
# pol vs. tor here except for the output plot. The resulting output is identical.
if pol_lims:
    wall_path = "/Users/zamperini/Documents/d3d_work/lwt/167196/mafot_3d_wall.dat"
    if shot == 167196:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/167196/167196_3500.pickle"
        output_dict = "ncoords_167196_pol.pickle"
        output_file = "for_mafot_167196_pol.dat"
    elif shot == 186914:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/186914_3500.pickle"
        output_dict = "ncoords_186914_pol.pickle"
        output_file = "for_mafot_186914_pol.dat"
    elif shot == 176971:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/176971_3000.pickle"
        output_dict = "ncoords_176971_pol.pickle"
        output_file = "for_mafot_176971_pol.dat"
else:
    wall_path = "/Users/zamperini/Documents/d3d_work/lwt/930116/mafot_wall_wide.dat"
    if shot == 167196:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/167196/167196_3500.pickle"
        output_dict = "ncoords_167196_tor.pickle"
        output_file = "for_mafot_167196_tor.dat"
    elif shot == 186914:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/186914_3500.pickle"
        output_dict = "ncoords_186914_tor.pickle"
        output_file = "for_mafot_186914_tor.dat"
    elif shot == 176971:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/176971_3000.pickle"
        output_dict = "ncoords_176971_tor.pickle"
        output_file = "for_mafot_176971_tor.dat"


# Load gfile.
with open(gfile_pickle_path, "rb") as f:
    gfile = pickle.load(f)
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
psin = gfile["PSIRZ_NORM"]
R_sep = gfile["RBBBS"]
Z_sep = gfile["ZBBBS"]

# Limit to just the right of the OMP.
raxis = gfile["RMAXIS"]
outboard = gfile["R"] > raxis
gR = gR[:, outboard]
gZ = gZ[:, outboard]
psin = psin[:, outboard]

# Grab wall coordinates.
lwt = LimWallToolkit()
wall = lwt.read_3d_wall(wall_path)

# Plot to show ranges of psin selected.
fig, ax = plt.subplots()
levels = [psins[0], psins[-1]]
ax.contour(gR, gZ, psin, levels=levels, colors="k", linewidths=1)
ax.contour(gR, gZ, psin, levels=[1.0], linewidths=3, colors="k")
ax.plot(wall[242][0], wall[242][1], lw=2, color="k")
ax.set_aspect("equal")
fig.tight_layout()
fig.show()

# Take advantage of contour to get the (R, Z) coordinates for each psin. A
# dictionary will be pickled that contains the number of coordinates per
# psin so that you can eventually back out of the MAFOT results which
# coordinates apply to which psin.
conts = ax.contour(gR, gZ, psin, levels=psins)
ncoords = {}
with open(output_file, "w") as f:
    f.write("# R, degree, Z\n")
    for p, seg in zip(psins, conts.allsegs):
        if len(seg) > 1:
            print("Error! More than one contour: {:.4f}".format(p))
            continue
        else:
            print("{:.4f}: {}".format(p, len(seg[0])))
        for deg in range(0, 360, 2):
            n = 0
            for i in range(0, len(seg[0])):

                # Only write every 5th data point on the segment.
                if i % 5 == 0:
                    f.write("{:.4f}\t{:3}\t{:.4f}\n".format(seg[0][i,0], deg, seg[0][i,1]))
                    n += 1
            if deg == 0:
                ncoords[p] = n

with open(output_dict, "wb") as f:
    pickle.dump(ncoords, f)
