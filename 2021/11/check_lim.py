# Script to check the output of a 3DLIM run to make sure it seems reasonable.
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import matplotlib as mpl


# Only options to change are here.
#ncpath = "/Users/zamperini/Documents/lim_runs/mcp4-184527-tor240_flatdif.nc"
#ncpath = "/Users/zamperini/Documents/d3d_work/167196/167196-a2-tor240_50a.nc"
#ncpath = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-005.nc"
#ncpath = "/Users/zamperini/Documents/d3d_work/184527/mcp4-184527-tor240_m0_v3.nc"
ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/tor_lim_testing/167196-pol-bump-002.nc"
pol_idx = 21  # 41 poloidal bins, 20 is the middle.
#plot_data = "velplasma2"
#plot_data = "velplasma_4d_1"  # Velocity in regions without CP.
#plot_data = "velplasma_4d_2"  # Velocity in regions with CP.
#plot_data = "crnbs_3d"        # Background plasma density.
#plot_data = "ctembs_3d"       # Background plasma temperature.
plot_data = "ddlim3"          # Impurity density.
ddlim3_charge = None

# Load netcdf file.
nc = netCDF4.Dataset(ncpath)

# Copied from lim_plots.
cl = float(nc['CL'][:].data)

# Same with the location of the plasma center (the top of the box)
ca  = float(nc['CA'][:].data)
caw = float(nc['CAW'][:].data)

# Get the X and Y grid data.
x = nc.variables['XOUTS'][:].data
y = nc.variables['YOUTS'][:].data
p = nc.variables["PS"][:].data
ca = nc.variables["CA"][:].data

# 2D representation of the data in the radial, parallel plane.
if plot_data == "ddlim3":
    if type(ddlim3_charge) == type(None):
        Z = nc.variables[plot_data.upper()][:].data[pol_idx, :, :, :].sum(axis=0)
    else:
        Z = nc.variables[plot_data.upper()][:].data[pol_idx, ddlim3_charge, :, :]
    Z2 = nc.variables[plot_data.upper()][:].data[:, :, :, :].sum(axis=1)
elif plot_data == "velplasma1":
    Z = nc.variables["velplasma"][:].data[0, :, :]
    Z2 = nc.variables["velplasma"][:].data
elif plot_data == "velplasma2":
    Z = nc.variables["velplasma"][:].data[1, :, :]
    Z2 = nc.variables["velplasma"][:].data
else:
    Z = nc.variables[plot_data][:].data[:, :, pol_idx]
    Z2 = nc.variables[plot_data][:].data[:, :, :]

# Trim the zeros from the edges of the x and y arrays, and the associated
# data points as well. This is done to stop this data from messing up
# the contours in the contour plot.
xkeep_min = np.nonzero(x)[0].min()
xkeep_max = np.nonzero(x)[0].max()
ykeep_min = np.nonzero(y)[0].min()
ykeep_max = np.nonzero(y)[0].max()
x = x[xkeep_min:xkeep_max]
y = y[ykeep_min:ykeep_max]
Z = Z[ykeep_min:ykeep_max, xkeep_min:xkeep_max]
if plot_data == "ddlim3":
    Z2 = Z2[:, ykeep_min:ykeep_max, xkeep_min:xkeep_max]
elif plot_data in["velplasma1", "velplasma2"]:
    Z2 = None
else:
    Z2 = Z2[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

# Furthermore, trim the data off that is beyond CL.
ykeep_cl = np.where(np.abs(y) < cl)[0]
y = y[ykeep_cl]
Z = Z[ykeep_cl, :]
if plot_data == "ddlim3":
    Z2 = Z2[:, ykeep_cl, :]
elif plot_data in["velplasma1", "velplasma2"]:
    Z2 = None
else:
    Z2 = Z2[ykeep_cl, :]

# Load in the bounds at this poloidal index.
try:
    bounds1a = nc.variables["bounds_1a"][:].data
    bounds2a = nc.variables["bounds_2a"][:].data
    bounds1a_slice = bounds1a[pol_idx]
    bounds2a_slice = bounds2a[pol_idx]
except:
    pass

# Mask zeros.
if plot_data in ["crnbs_3d", "ctembs_3d"]:
    Z = np.ma.masked_where(Z<=0, Z)
    Z2 = np.ma.masked_where(Z2<=0, Z2)
    cmap = "inferno"
    vmin = Z.min()
    vmax = Z.max()
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
elif plot_data in ["velplasma1", "velplasma2", "velplasma_4d_1", "velplasma_4d_2"]:
    cmap = "coolwarm"
    vmin = -np.abs(Z).max()
    vmax = np.abs(Z).max()
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
elif plot_data == "ddlim3":
    Z = np.ma.masked_where(Z<=0, Z)
    Z2 = np.ma.masked_where(Z2<=0, Z2)
    cmap = "inferno"
    vmin = Z.max() / 1000
    vmax = Z.max() / 1
    norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)

Y, X = np.meshgrid(x, y)

fig, ax = plt.subplots()

#cont = ax.contourf(X, Y, Z, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm)
cont = ax.pcolormesh(X, Y, Z, shading="auto", norm=norm, cmap=cmap)
cbar = fig.colorbar(cont)
#try:
ax.step(bounds1a_slice, np.append(x, ca), color="k", where="post")
ax.step(bounds2a_slice, np.append(x, ca), color="k", where="post")
#except:
#    pass
ax.set_xlabel("Parallel to B (m)", fontsize=14)
ax.set_ylabel("Radial (m)", fontsize=14)
cbar.set_label("C Density (arbitrary)", fontsize=14)

fig.tight_layout()
fig.show()

# Top-down ddlims3.
#if plot_data not in ["velplasma1", "velplasma2"]:
#    fig, ax = plt.subplots()
#    if plot_data == "ddlim3":
#        Z2_plot = Z2[:, :, -3]
#    else:
#        Z2_plot = Z2[:, -3, :].T
#    cont = ax.pcolormesh(y, p, Z2_plot, shading="auto", vmin=Z2_plot.max()/10, vmax=Z2_plot.max(), cmap=cmap)
#    ax.set_ylabel("P (m)")
#    ax.set_xlabel("Y (m)")
#    fig.tight_layout()
#    fig.show()

# Plot of the density along the top few rows.
fig, ax = plt.subplots()
ax.plot(y, Z[:,-1])
ax.plot(y, Z[:,-2])
ax.plot(y, Z[:,20])
ax.plot(y, Z[:,30])
ax.set_xlabel("Y (m)")
ax.set_ylabel("DDLIM3")
fig.tight_layout()
fig.show()

# Show the bounds.
try:
    vmax = np.max((np.abs(bounds1a), np.abs(bounds2a)))
    fig, (ax1, ax2) = plt.subplots(1, 2)
    cont1 = ax1.pcolormesh(p, x, bounds1a.T, shading="auto", vmax=vmax)
    cont2 = ax2.pcolormesh(p, x, bounds2a.T*-1, shading="auto", vmax=vmax)
    ax1.axvline(p[pol_idx], color="r", lw=2)
    ax2.axvline(p[pol_idx], color="r", lw=2)
    ax1.set_title("Positive Bounds")
    ax2.set_title("Negative Bounds")
    cbar = fig.colorbar(cont2, ax=ax2)
    fig.tight_layout()
    fig.show()
except:
    pass
