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
#psins = np.linspace(1.18, 1.45, 50)   # 186914
#psins = np.linspace(1.20, 1.50, 50)   # 167463
shot = 190422
pol_lims = True  # Leave True

# Choose correct paths for everything. There actually is no reason to separate
# pol vs. tor here except for the output plot. The resulting output is identical.
if pol_lims:
    wall_path = "/Users/zamperini/Documents/d3d_work/lwt/167196/mafot_3d_wall.dat"
    if shot == 167196:
        psins = np.linspace(1.13, 1.40, 50)
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/ncoords_167196_pol.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/for_mafot_167196_pol.dat"
    elif shot == 167463:
        psins = np.linspace(1.20, 1.50, 50)
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/167463_3000.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/ncoords_167463_pol.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/for_mafot_167463_pol.dat"
    elif shot == 180455:
        psins = np.linspace(1.20, 1.50, 50)
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/180455/180455_3500.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot_files/180455/ncoords_180455.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot_files/180455/for_mafot_180455.dat"
    elif shot == 186257:
        psins = np.linspace(1.20, 1.45, 50)
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/186257/186257_3500.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot_files/186257/ncoords_186257.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot_files/186257/for_mafot_186257.dat"
    elif shot == 186754:
        psins = np.linspace(1.18, 1.45, 50)
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/186754/186754_3500.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot_files/186754/ncoords_186754.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot_files/186754/for_mafot_186754.dat"
    elif shot == 186914:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/186914_3500.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/ncoords_186914_pol.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/for_mafot_186914_pol.dat"
    elif shot == 176971:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/176971_3000.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/ncoords_176971_pol.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/for_mafot_176971_pol.dat"
    elif shot == 174783:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/174783_0.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/ncoords_174783_svr.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/for_mafot_174783_svr.dat"
    elif shot == 190422:
        psins = np.linspace(1.20, 1.50, 50)
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190422/190422_3500.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot_files/190422/ncoords_190422.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot_files/190422/for_mafot_190422.dat"
else:
    wall_path = "/Users/zamperini/Documents/d3d_work/lwt/930116/mafot_wall_wide.dat"
    if shot == 167196:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot/167196/167196_3500.pickle"
        output_dict = "/Users/zamperini/Documents/d3d_work/mafot/167196/ncoords_167196_tor.pickle"
        output_file = "/Users/zamperini/Documents/d3d_work/mafot/167196/for_mafot_167196_tor.dat"
    elif shot == 167463:
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/167463_3000.pickle"
        output_dict = "ncoords_167463_tor.pickle"
        output_file = "for_mafot_167463_tor.dat"
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
