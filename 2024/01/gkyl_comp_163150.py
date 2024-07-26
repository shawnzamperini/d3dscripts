import postgkyl as pg
import matplotlib.pyplot as plt
import numpy as np
import pickle
from freeqdsk import geqdsk


# Geometry information from input file.
Rdim = 1.7  # [m]
Zdim = 3.2  # [m]
Z_axis = 0.00232616113
R_axisTrue = 1.72068012
R_axis = R_axisTrue
B_axis = 2 * R_axisTrue / R_axis
R_LCFSmid = 2.2801477223421736
Rmid_min = R_LCFSmid - 0.1
Rmid_max = R_LCFSmid + 0.05
R0 = 0.5 * (Rmid_min + Rmid_max)
a_shift = 0.6
a_mid = R_axis / a_shift - np.sqrt(R_axis * (R_axis - 2 * a_shift * R_LCFSmid + 2 * a_shift * R_axis)) / a_shift
r0 = R0 - R_axis
B0 = B_axis * (R_axis / R0)
kappa = 1.488  # Elongation.
delta = 0.1  # Triangularity.

# Will want to map to the outboard midplane.
gfile_path = "/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-163150-v1/g163150.01500"
with open(gfile_path, "r") as f:
    gfile = geqdsk.read(f)
gfile_R = np.linspace(gfile["rleft"], gfile["rleft"] + gfile["rdim"], gfile["nx"])
gfile_Z = np.linspace(gfile["zmid"] - gfile["zdim"] / 2, gfile["zmid"] + gfile["zdim"] / 2, gfile["ny"])
gfile_Rs, gfile_Zs = np.meshgrid(gfile_R, gfile_Z, indexing="ij")
gfile_psis = gfile["psi"]
rlim = gfile["rlim"]
zlim = gfile["zlim"]

# Functions to map gkyl coordinates to (R, Z).

def r_x(xIn):
    """
    Minor radius as a function of x.
    """
    return Rmid_min+xIn-R_axis


def R_f(r, theta):
    return R_axis + r * np.cos(theta + np.arcsin(delta) * np.sin(theta))


def Z_f(r, theta):
    return kappa * r * np.sin(theta)


# Adapted from the postgkyl documentation: https://gkeyll.readthedocs.io/en/latest/postgkyl/examplesScript.html
def load_gkyl_data_rad(path, z_slice):
    data = pg.data.GData(path)
    data_interp = pg.data.GInterpModal(data, 1, 'ms')
    interp_grid, interp_values = data_interp.interpolate()

    # get cell center coordinates
    CCC = []
    for j in range(0, len(interp_grid)):
        CCC.append((interp_grid[j][1:] + interp_grid[j][:-1]) / 2)

    # Axisymmetric, so no y values.
    x_vals = CCC[0]
    z_vals = CCC[1]

    print("len(z_vals) = {}".format(len(z_vals)))
    print("z_slice at {:.2f} degrees".format(z_vals[z_slice] * 180 / np.pi))

    X, Z = np.meshgrid(x_vals, z_vals)
    rad_data = np.transpose(interp_values[:, z_slice, 0])
    return x_vals, rad_data


def load_gkyl_data_along_line(path, r_start, z_start, r_end, z_end, n_line_points=100, show_plot=False, scale=1.0,
                              comp=0):
    """
    Load data from gkyl using data along a line that goes from (r_start, z_start) to (r_end, z_end).
    """
    data = pg.data.GData(path)
    data_interp = pg.data.GInterpModal(data, 1, 'ms')
    interp_grid, interp_values = data_interp.interpolate(comp)
    # gkyl_data = np.transpose(interp_values[:, :, 0])
    gkyl_data = interp_values[:, :, 0]
    gkyl_data *= scale

    # get cell center coordinates
    CCC = []
    for j in range(0, len(interp_grid)):
        CCC.append((interp_grid[j][1:] + interp_grid[j][:-1]) / 2)

    # Axisymmetric, so no y values. z in this instance is theta, not Z in machine coordinates. That happens next.
    gkyl_x_vals = CCC[0]
    gkyl_theta_vals = CCC[1]
    gkyl_X, gkyl_Theta = np.meshgrid(gkyl_x_vals, gkyl_theta_vals)

    # Convert to (R, Z) coordinates.
    gkyl_R = np.zeros((len(gkyl_x_vals), len(gkyl_theta_vals)))
    gkyl_Z = np.zeros(gkyl_R.shape)
    for i in range(gkyl_R.shape[0]):
        for j in range(gkyl_R.shape[1]):
            gkyl_R[i, j] = R_f(r_x(gkyl_x_vals[i]), gkyl_theta_vals[j])
            gkyl_Z[i, j] = Z_f(r_x(gkyl_x_vals[i]), gkyl_theta_vals[j])

    if show_plot:
        fig, ax = plt.subplots()
        ax.set_aspect("equal")
        ax.plot(rlim, zlim, color="k")
        ax.contourf(gkyl_R, gkyl_Z, gkyl_data)
        ax.plot([r_start, r_end], [z_start, z_end], color="k", lw=4)
        ax.plot([r_start, r_end], [z_start, z_end], color="r", lw=3)
        ax.scatter(gkyl_R, gkyl_Z, color="w", edgecolors="k", s=3)
        fig.tight_layout()
        fig.show()

    # Point-slope form for the line we are pulling data along.
    line_rs = np.linspace(r_start, r_end, n_line_points)
    if r_start != r_end:
        line_slope = (z_end - z_start) / (r_end - r_start)
        line_zs = line_slope * (line_rs - r_start) + z_start
    else:
        line_zs = np.linspace(z_start, z_end, n_line_points)

    # For each point along the line, pull the closest gkyl data point and store it.
    line_stored_r = []
    line_stored_z = []
    line_stored_data = []
    stored_idx = []
    for n in range(n_line_points):
        dist = np.sqrt(np.square(line_rs[n] - gkyl_R) + np.square(line_zs[n] - gkyl_Z))
        min_idx = np.where(dist == dist.min())
        # print(min_idx)

        # Avoid storing duplicates.
        if min_idx in stored_idx:
            continue

        stored_idx.append(min_idx)
        line_stored_r.append(line_rs[n])
        line_stored_z.append(line_zs[n])
        line_stored_data.append(gkyl_data[min_idx][0])


    return line_stored_r, line_stored_z, line_stored_data


# Load TS data.
gfile_path = "/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-163150-v1/ts_163150.pickle"
with open(gfile_path, "rb") as f:
    ts = pickle.load(f)

time = 1500
trange = 100
# rsep = 1.94
# zsep = 0.473
# idx = np.argmin(np.abs(ts["core"]["time"] - time))
idx = np.logical_and(ts["core"]["time"] > (time - trange), ts["core"]["time"] < time + trange)
psin = ts["core"]["psin"][:, idx].flatten()
te = ts["core"]["te"][:, idx].flatten()
ne = ts["core"]["ne"][:, idx].flatten()
# r = np.tile(ts["core"]["r"], idx.sum())
# z = np.tile(ts["core"]["z"], idx.sum())
r = np.repeat(ts["core"]["r"], idx.sum())
z = np.repeat(ts["core"]["z"], idx.sum())

# Not interested in anything where psin = 0.
mask = psin != 0
psin = psin[mask]
te = te[mask]
ne = ne[mask]
r = r[mask]
z = z[mask]

# Sort by psin value.
sort_idx = np.argsort(psin)
psin = psin[sort_idx]
te = te[sort_idx]
ne = ne[sort_idx]
r = r[sort_idx]
z = z[sort_idx]

# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# ax.contourf(gfile_Rs, gfile_Zs, gfile_psis)
# ax.contour(gfile_Rs, gfile_Zs, gfile_psis, levels=[gfile["sibdry"]], colors="r")
# ax.plot(rlim, zlim, color="k")
# fig.tight_layout()
# fig.show()

# Now load the CER data that looks at C6+.
cer_path = "/Users/zamperini/github/d3dscripts/2024/01/cer_163150_1400_1600.pickle"
with open(cer_path, "rb") as f:
    cer = pickle.load(f)

# Put into convenient arrays.
cer_r = np.array([cer[k]["r"] for k in cer.keys()])
cer_nz = np.array([cer[k]["nz"] for k in cer.keys()])
cer_tz = np.array([cer[k]["tz"] for k in cer.keys()])

# Load BES data. Note this doesn't matter at this point because axisymmetric = no fluctuations.

# TS crosses separatrix at (R, Z) = (1.94 0.792). Magnetic axis at (1.721, 0.002). That's:
#                                / |
#                          /       |
#                     /            | 0.79 m
#                /                 |
#           / a                    |
#        ---------------------------
#             0.219 m
#
#  tan(a) = 0.79 / 0.219, or a = 1.30 rad = 74.5 degrees.
elc_ne_path = "/Users/zamperini/gkyldir/d3d-163150-axisym-v2/d3d-163150-axisym-v2-elc_M0_100.gkyl"
elc_te_path = "/Users/zamperini/gkyldir/d3d-163150-axisym-v2/d3d-163150-axisym-v2-elc_prim_moms_100.gkyl"
imp_nz_path = "/Users/zamperini/gkyldir/d3d-163150-axisym-v2/d3d-163150-axisym-v2-imp_M0_100.gkyl"
imp_tz_path = "/Users/zamperini/gkyldir/d3d-163150-axisym-v2/d3d-163150-axisym-v2-imp_prim_moms_100.gkyl"
gkyl_ne_z0, gkyl_ne_rad = load_gkyl_data_rad(elc_ne_path, z_slice=45)  # z_slice=45 is about 76 degrees
gkyl_nz_z0, gkyl_nz_rad = load_gkyl_data_rad(imp_nz_path, z_slice=0)   # CER is at the outboard midplane, 0 degrees

# Get gkyl data along span of TS.
r_start = 1.94  # Actual line that TS covers.
z_start = 0.56
r_end = 1.94
z_end = 0.825
gkyl_ne_r, gkyl_ne_z, gkyl_ne = load_gkyl_data_along_line(elc_ne_path, r_start=r_start, z_start=z_start, r_end=r_end,
                                                          z_end=z_end, show_plot=True)
gkyl_te_r, gkyl_te_z, gkyl_te = load_gkyl_data_along_line(elc_te_path, r_start=r_start, z_start=z_start, r_end=r_end,
                                                          z_end=z_end, show_plot=False, scale=9.1e-31/1.602e-19,
                                                          comp=1)

# gkyl C6+ data at the outboard midplane
r_start = 2.16
z_start = 0.01  # Small offset because Z=0 is equidistant between two cells.
r_end = 2.32
z_end = 0.01
gkyl_nz_r, gkyl_nz_z, gkyl_nz = load_gkyl_data_along_line(imp_nz_path, r_start=r_start, z_start=z_start, r_end=r_end,
                                                          z_end=z_end, show_plot=True)

# Before doing carbon, need the carbon mass so we can go from vTz --> Tz.
mc = 12.011 * 1.66e-27  # kg
gkyl_tz_r, gkyl_tz_z, gkyl_tz = load_gkyl_data_along_line(imp_tz_path, r_start=r_start, z_start=z_start, r_end=r_end,
                                                          z_end=z_end, show_plot=False, scale=mc/1.602e-19,
                                                          comp=1)

# Te is actually the thermal velocity, so convert to Te. 1.602e-19 / 9.1e-31 *
gkyl_te = np.array(gkyl_te)

# Plot comparing Gkeyll to TS data.
zsep = 0.794
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 7))

# Density
ax1.axvline(zsep, color="k")
ax1.scatter(z, ne, edgecolors="k", color="tab:red")
ax1.plot(gkyl_ne_z, gkyl_ne, color="k", lw=4)
ax1.plot(gkyl_ne_z, gkyl_ne, color="tab:red", lw=3)
ax1.set_xlabel("Z (m)", fontsize=12)
ax1.set_ylabel("ne (m-3)", fontsize=12)
ax1.set_xlim([0.5, None])

# Temperature
ax2.axvline(zsep, color="k")
ax2.scatter(z, te, edgecolors="k", color="tab:red")
ax2.plot(gkyl_te_z, gkyl_te, color="k", lw=4)
ax2.plot(gkyl_te_z, gkyl_te, color="tab:red", lw=3)
ax2.set_xlabel("Z (m)", fontsize=12)
ax2.set_ylabel("Te (eV)", fontsize=12)
ax2.set_ylim([0, 500])
ax2.set_xlim([0.5, None])

# C6+ density.
ax3.axvline(2.28, color="k")
ax3.scatter(cer_r, cer_nz, edgecolors="k", color="tab:green")
ax3.plot(gkyl_nz_r, gkyl_nz, color="k", lw=4)
ax3.plot(gkyl_nz_r, gkyl_nz, color="tab:green", lw=3)
ax3.set_xlabel("R (m)", fontsize=12)
ax3.set_ylabel("C6+ Density (m-3)", fontsize=12)
ax3.set_xlim([2.1, None])

# C6+ temperature.
ax4.axvline(2.28, color="k")
ax4.scatter(cer_r, cer_tz, edgecolors="k", color="tab:green")
ax4.plot(gkyl_tz_r, gkyl_tz, color="k", lw=4)
ax4.plot(gkyl_tz_r, gkyl_tz, color="tab:green", lw=3)
ax4.set_xlabel("R (m)", fontsize=12)
ax4.set_ylabel("C6+ Temperature (eV)", fontsize=12)
ax4.set_xlim([2.1, None])
ax4.set_ylim([0, 1000])

fig.tight_layout()
fig.show()