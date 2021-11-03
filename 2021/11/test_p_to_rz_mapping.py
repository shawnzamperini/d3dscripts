# Meant to be a test script to see how to go from specified 3DLIM bins to an
# R, Z coordinate required for getting the connection lengths.
import pickle
import numpy as np
import matplotlib.pyplot as plt


# 3DLIM R and P bins.
lim_rbins = np.arange(-0.1, 0.02, 0.0025)
lim_pbins = np.linspace(-0.8, 0.8, 41)

mimes = False
if mimes:
    # MiMES style origin.
    r_origin2 = 2.295
    z_origin2 = -0.188
else:
    # DiMES style origin.
    r_origin2 = 1.485
    z_origin2 = -1.205  # 4 cm above the floor.

# gfile with equilibrium.
gfile_pickle_path2 = "/Users/zamperini/Documents/d3d_work/184527/184527_3500.pickle"
with open(gfile_pickle_path2, "rb") as f:
    gfile = pickle.load(f)
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
psin  = gfile["PSIRZ_NORM"]
R_sep = gfile["RBBBS"]
Z_sep = gfile["ZBBBS"]


def get_psin(r, z):
    dist = np.sqrt(np.square(r-gR) + np.square(z-gZ))
    close_psin = psin[np.where(dist == dist.min())]

    return {"close_psin":close_psin[0]}

def get_mach_coords(mid_r, mid_z, lim_pbins):

    dist = np.sqrt(np.square(mid_r-gR) + np.square(mid_z-gZ))
    psin_mid = psin[np.where(dist == dist.min())]

    fig, ax = plt.subplots()

    # Include the correct half of the flux tubes (i.e. avoid outside the vessel).
    if mimes:
        keep = gR > 1.50
        keep = keep[0,:]
        gR_keep = gR[:,keep]
        gZ_keep = gZ[:,keep]
        psin_keep = psin[:,keep]
    else:
        #if usn:
        #    keep = gZ < 1.00
        #    keep = keep[:,0]
        #else:
        #    keep = gZ > 1.35
        #    keep = keep[:,0]
        keep = gZ < 0.00
        keep = keep[:,0]
        gR_keep = gR[keep]
        gZ_keep = gZ[keep]
        psin_keep = psin[keep]

    cont2 = ax.contour(gR_keep, gZ_keep, psin_keep, levels=[psin_mid], colors="r")
    mid_line_rz = cont2.allsegs[0][0]

    dist = np.sqrt(np.square(mid_line_rz[:,0] - mid_r) + np.square(mid_line_rz[:,1] - mid_z))
    close_mid_idx = np.where(dist == dist.min())[0][0]
    close_mid = mid_line_rz[close_mid_idx]

    # To hold the corresponding R, Z coordinate for each P bin.
    lim_machR = np.zeros(lim_pbins.shape)
    lim_machZ = np.zeros(lim_pbins.shape)

    pbin_center_idx = np.where(lim_pbins == 0)[0][0]

    # Go in one direction from the origin along the field line.
    pbin_idx = pbin_center_idx
    tot_dist = 0
    prev_coord = close_mid
    for d in mid_line_rz[close_mid_idx:]:

        # Tally distance travelled along field line, assigning an R, Z coordinate
        # each time we pass a P bin.
        tmp_dist = np.sqrt(np.square(d[0] - prev_coord[0]) + np.square(d[1] - prev_coord[1]))
        tot_dist += tmp_dist
        if tot_dist >= np.abs(lim_pbins[pbin_idx]):
            lim_machR[pbin_idx] = d[0]
            lim_machZ[pbin_idx] = d[1]
            pbin_idx += 1
            if pbin_idx > len(lim_pbins) - 1:
                break
        prev_coord = d

    # Likewise in the other direction.
    pbin_idx = pbin_center_idx
    tot_dist = 0
    prev_coord = close_mid
    for d in mid_line_rz[:close_mid_idx][::-1]:
        tmp_dist = np.sqrt(np.square(d[0] - prev_coord[0]) + np.square(d[1] - prev_coord[1]))
        tot_dist += tmp_dist
        if tot_dist >= np.abs(lim_pbins[pbin_idx]):
            lim_machR[pbin_idx] = d[0]
            lim_machZ[pbin_idx] = d[1]
            pbin_idx -= 1
            if pbin_idx < 0:
                break
        prev_coord = d

    plt.close()
    return {"lim_machR":lim_machR, "lim_machZ":lim_machZ}

fig, ax = plt.subplots()
cont1 = ax.contour(gR, gZ, psin, colors="k", levels=np.arange(0.95, 1.10, 0.025), linewidths=1)
cont1 = ax.contour(gR, gZ, psin, colors="k", levels=[1], linewidths=3)

for lim_r in lim_rbins:

    if mimes:
        # Get the corresponding machine R coordinate.
        machR = r_origin2 - lim_r
        machZ = z_origin2
    else:
        # If DiMES, then the 3DLIM R coord is actually parallel to Z.
        machZ = z_origin2 - lim_r
        machR = r_origin2

    mach_coords = get_mach_coords(machR, machZ, lim_pbins)
    lim_machR = mach_coords["lim_machR"]
    lim_machZ = mach_coords["lim_machZ"]
    ax.scatter(lim_machR, lim_machZ, color="r", s=10)
ax.set_aspect("equal")
fig.tight_layout()
fig.show()



"""
# Get the R, Z coordinates of this psin contour.
fig, ax = plt.subplots()
cont1 = ax.contour(gR, gZ, psin, colors="k")

# Need to do the below once for each flux tube we are on, represented by the
# R-Rsep coordinate of the lim_rbins.
if mimes:
    lim_rbins_machR = lim_rbins - r_origin2
else:
    lim_rbins_machR = np.full()

# Identify the psin of the flux surface we are on.
psin_origin = get_psin(r_origin2, z_origin2)["close_psin"]
cont2 = ax.contour(gR, gZ, psin, levels=[psin_origin], colors="r")

origin_line_rz = cont2.allsegs[0][0]

# Identify where on our flux line is closest to the origin.
dist = np.sqrt(np.square(origin_line_rz[:,0] - r_origin2) + np.square(origin_line_rz[:,1] - z_origin2))
close_origin_idx = np.where(dist == dist.min())[0][0]
close_origin = origin_line_rz[close_origin_idx]

# To hold the corresponding R, Z coordinate for each P bin.
lim_machR = np.zeros(lim_pbins.shape)
lim_machZ = np.zeros(lim_pbins.shape)

pbin_center_idx = np.where(lim_pbins == 0)[0][0]

# Go in one direction from the origin along the field line.
pbin_idx = pbin_center_idx
tot_dist = 0
prev_coord = close_origin
for d in origin_line_rz[close_origin_idx:]:

    # Tally distance travelled along field line, assigning an R, Z coordinate
    # each time we pass a P bin.
    tmp_dist = np.sqrt(np.square(d[0] - prev_coord[0]) + np.square(d[1] - prev_coord[1]))
    tot_dist += tmp_dist

    if tot_dist >= np.abs(lim_pbins[pbin_idx]):
        lim_machR[pbin_idx] = d[0]
        lim_machZ[pbin_idx] = d[1]
        pbin_idx += 1
        if pbin_idx > len(lim_pbins) - 1:
            break

    prev_coord = d

# Likewise in the other direction.
pbin_idx = pbin_center_idx
tot_dist = 0
prev_coord = close_origin
for d in origin_line_rz[:close_origin_idx][::-1]:

    # Tally distance travelled along field line, assigning an R, Z coordinate
    # each time we pass a P bin.
    tmp_dist = np.sqrt(np.square(d[0] - prev_coord[0]) + np.square(d[1] - prev_coord[1]))
    tot_dist += tmp_dist

    if tot_dist >= np.abs(lim_pbins[pbin_idx]):
        lim_machR[pbin_idx] = d[0]
        lim_machZ[pbin_idx] = d[1]
        pbin_idx -= 1
        if pbin_idx < 0:
            break

    prev_coord = d

ax.scatter(lim_machR, lim_machZ)
ax.set_aspect("equal")
fig.tight_layout()
fig.show()
"""
