# Script to plot the sputtered impurity density around a limiter in 3DLIM. Only
# half of the limiter is shown due to the design choices in 3DLIM. May be able
# to eventually plot both sides of the limiter one day.
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import matplotlib as mpl
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.lines import Line2D


#ncpath = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-010-otf.nc"
ncpath = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-5cm-lim-001.nc"
pol_idx = 21
#plot_data = "velplasma1"     # Velocity in regions without (1) or with (1) limiter.
#plot_data = "velplasma_4d_1"  # Velocity in regions without CP.
#plot_data = "velplasma_4d_2"  # Velocity in regions with CP.
#plot_data = "crnbs_3d"        # Background plasma density.
#plot_data = "ctembs_3d"       # Background plasma temperature.
#plot_data = "CTEMBS"
plot_data = "ddlim3"          # Impurity density.
ddlim3_charge = None             # 2 = C+, 3 = C2+, ..., 7 = C6+
include_neutrals = False

# Window for the plot so we can zoom in on the limiter.
xlim1 = [-0.2, 0.2]; equal = True
#xlim1 = [-45, 45]; equal = False

# Load netcdf file and some of the output data.
print("Loading in {}...".format(ncpath))
nc = netCDF4.Dataset(ncpath)
cl = float(nc['CL'][:].data)
ca  = float(nc['CA'][:].data)
caw = float(nc['CAW'][:].data)
xs = nc.variables["XS"][:].data
ys = nc.variables["YS"][:].data
ps = nc.variables["PS"][:].data
xwids = nc.variables["XWIDS"][:].data
ywids = nc.variables["YWIDS"][:].data
pwids = nc.variables["PWIDS"][:].data
xouts = nc.variables['XOUTS'][:].data
youts = nc.variables['YOUTS'][:].data
absfac = nc.variables["ABSFAC"][:].data
crnbs = nc.variables["CRNBS"][:].data
odouts = nc.variables['ODOUTS'][:].data
nerods3 = nc.variables["NERODS3"][:].data[1]  # 0 = deposition. 1 = primary erosion. 2 = erosion + self-sputtering. 4 = net erosion
nerods = nc.variables["NERODS"][:].data[1]
qtembs = nc.variables["QTEMBS"][:][0][:-2]  # Last 2 values seem to be junk.
odouts = nc.variables["ODOUTS"][:]
qtans = nc.variables["QTANS"][:]
qtembsi = nc.variables["QTEMBSI"][:][0][:-2]

print("ABSFAC = {:.2e}".format(absfac))

# NERODS3 has units of particles/m2/# of particles

# Need to calculate pouts.
pouts = ps - pwids / 2

# Mirror the Y data around 0.
ys = np.append(np.append(-ys[::-1], 0), ys)
ywids = np.append(np.append(-ywids[::-1], 0), ywids)
youts = np.append(np.append(-youts[::-1], 0), youts)

# The limiter locations.
qxs = nc.variables["QXS"][:].data
qedges1 = nc.variables["QEDGES"][:].data[0]
qedges2 = nc.variables["QEDGES"][:].data[1]

# The D ion impact energy.
eimp = 2 * qtembsi + 3 * qtembs

#qedges_tmp = nc.variables["QEDGES"][:].data
#fig, ax1 = plt.subplots()
#ax11 = ax1.twinx()
#ax1.plot(qxs, "k")
#ax11.plot(qedges_tmp[0], "r")
#ax11.plot(qedges_tmp[1], "b")
#fig.tight_layout()
#fig.show()

# 2D representation of the data in the radial, parallel plane.
if plot_data == "ddlim3":
    if type(ddlim3_charge) == type(None):

        # The first index of the second dimension is primary neutrals (0), while
        # the second is all neutrals (1). It then goes into the ionization states
        # after that (2+).
        if include_neutrals:
            Z = nc.variables[plot_data.upper()][:].data[pol_idx, 1:, :, :].sum(axis=0)
        else:
            Z = nc.variables[plot_data.upper()][:].data[pol_idx, 2:, :, :].sum(axis=0)
    else:
        Z = nc.variables[plot_data.upper()][:].data[pol_idx, ddlim3_charge, :, :]
    Z = Z * absfac
    Z = np.ma.masked_where(Z<=0, Z)
    norm = LogNorm(vmin=1e14, vmax=1e18)

elif plot_data in ["CTEMBS"]:
    Z = nc.variables[plot_data][:].data
    norm = None

else:
    Z = nc.variables[plot_data][:].data[:, :, pol_idx]
    norm = None

# If bounds_1a is in here then vary_2d_bounds was on.
if "bounds_1a" in nc.variables.keys():
    vp1 = nc.variables["velplasma_4d_1"][:, :, pol_idx]
    vp2 = nc.variables["velplasma_4d_2"][:, :, pol_idx]
else:
    vp1 = nc.variables["velplasma"][:].data[0, :, :]
    vp2 = nc.variables["velplasma"][:].data[1, :, :]

# We then need to trim all the trailing zeros on the arrays.
xkeep_min = np.nonzero(xs)[0].min()
xkeep_max = np.nonzero(xs)[0].max()
ykeep_min = np.nonzero(ys)[0].min()
ykeep_max = np.nonzero(ys)[0].max()
xs = xs[xkeep_min:xkeep_max]
ys = ys[ykeep_min:ykeep_max]
xwids = xwids[xkeep_min:xkeep_max]
ywids = ywids[ykeep_min:ykeep_max]
xouts = xouts[xkeep_min:xkeep_max]
youts = youts[ykeep_min:ykeep_max]
Z = Z[ykeep_min:ykeep_max, xkeep_min:xkeep_max]
crnbs = crnbs[ykeep_min:ykeep_max, xkeep_min:xkeep_max]
vp1 = vp1[ykeep_min:ykeep_max, xkeep_min:xkeep_max]
vp2 = vp2[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

#qxskeep_min = np.nonzero(qxs)[0].min()
#qxskeep_max = np.nonzero(qxs)[0].max()
#qxs = qxs[qxskeep_min:qxskeep_max]
#qedges = qedges[qxskeep_min:qxskeep_max]

# Only need the values up to len(qedges). Trim zeros. Want one zero at the
# end for fill between.
nonzero = np.nonzero(qedges1)
qxs = qxs[:len(qedges1)][nonzero]
qedges1 = qedges1[nonzero]
qedges2 = qedges2[nonzero]
qxs = np.append(qxs, 0)
qedges1 = np.append(qedges1, 0)
qedges2 = np.append(qedges2, 0)

# Set up for plotting.
half = int(len(qtembs)/2)
qtembs_plot = np.append(qtembs[half:], qtembs[half:][::-1])
qtans_plot = np.append(qtans[0][:-1][half:], qtans[1][half:][:-1][::-1])
eimp_plot = np.append(eimp[half:], eimp[half:][::-1])

# For the X, Y 2D arrays for pcolormesh, the documentation says these need to
# define the corner of each rectangle in the mesh, and should be one larger in
# each dimension. So go ahead and add on the last corners.
#x_coords = np.append(xs, xs[-1] + xwids[-1])
#y_coords = np.append(ys, ys[-1] + ywids[-1])
#x_coords = np.append(xs[0] + xwids[0], xs)
#y_coords = np.append(ys[0] + ywids[0], ys)
x_coords = np.append(xouts, xouts[-1] + xwids[-1] / 2)
y_coords = np.append(youts, youts[-1] + ywids[-1] / 2)
X, Y = np.meshgrid(x_coords, y_coords)

# 2D data for the nerods3 plot.
O, P = np.meshgrid(odouts, pouts)

# Index where X coordinate = 1 cm above the baffle and the indexed data.
#above_idx = np.argwhere(x_coords == 0.01)[0][0]
above_idx = 30
above_ddlim = Z[:, above_idx]

# Mask the nerods3 data.
nerods3_masked = np.ma.masked_where(nerods3<=0, nerods3)

# Set up plot grid of a respectable size.
fig, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(2, 3, figsize=(13, 7))

# Plot of our chosen data plotted around the limiter on a log scale.
cmesh = ax1.pcolormesh(Y, X, Z, cmap="inferno", norm=norm)
ax1.fill_between(-qedges1, qxs, X.min(), color="grey")
ax1.fill_between(qedges2, qxs, X.min(), color="grey")
ax1.plot(-qedges1, qxs, color="grey", lw=2)
ax1.plot(qedges2, qxs, color="grey", lw=2)
ax1.set_xlim(xlim1)
ax1.set_ylim([X.min(), ca])
if equal:
    ax1.set_aspect("equal")
cbar = fig.colorbar(cmesh, ax=ax1)
cbar.set_label("C Density (m-3)")
ax1.set_xlabel("Parallel to B (m)")
ax1.set_ylabel("Radial (m)")

# Plot of the impurity density as a fraction of the density.
cmesh = ax2.pcolormesh(Y, X, Z/crnbs, cmap="inferno", norm=LogNorm(vmin=0.001, vmax=1))
ax2.fill_between(-qedges1, qxs, X.min(), color="grey")
ax2.fill_between(qedges2, qxs, X.min(), color="grey")
ax2.plot(-qedges1, qxs, color="grey", lw=2)
ax2.plot(qedges2, qxs, color="grey", lw=2)
ax2.set_xlim(xlim1)
ax2.set_ylim([X.min(), ca])
if equal:
    ax2.set_aspect("equal")
cbar = fig.colorbar(cmesh, ax=ax2)
cbar.set_label("C Concentration")
ax2.set_xlabel("Parallel to B (m)")
ax2.set_ylabel("Radial (m)")

# Plot of erosion along the limiter with Te. Only index the negative valued
# surface until we can figure out what we're looking at.
side = int(len(odouts)/2)
ax33 = ax3.twinx()
#ax3.plot(-odouts[:side], -nerods[:side], label="Erosion", color="k")
ax3.plot(odouts, nerods3[pol_idx]*absfac, label="Erosion", color="k")
#ax33.plot(odouts, qtembs_plot, label="Te", color="tab:red")
ax33.plot(odouts, eimp_plot, label="Eimpact", color="tab:red")
#ax33.plot(odouts, qtans_plot, label="Te", color="tab:red")
custom_lines = [Line2D([0], [0], color="k", lw=2),
  Line2D([0], [0], color="tab:red", lw=2)]
ax3.legend(custom_lines, ["Erosion", "Te"], loc="upper right")
ax3.set_xlabel("Distance along limiter (m)")
ax3.set_ylabel("Net erosion (C/m2/s)")
ax33.set_ylabel(r"$\mathdefault{E_{impact}}$ at limiter (eV)", color="tab:red")
#ax3.set_xlim([0, -caw])
ax33.tick_params(axis='y', labelcolor="tab:red")

# Additional plot with a plasma background parameter of choice.
cmesh = ax4.pcolormesh(Y, X, vp2, cmap="coolwarm", norm=None)
ax4.fill_between(-qedges1, qxs, X.min(), color="grey")
ax4.fill_between(qedges2, qxs, X.min(), color="grey")
ax4.plot(-qedges1, qxs, color="grey", lw=2)
ax4.plot(qedges2, qxs, color="grey", lw=2)
ax4.set_xlim([-cl, cl])
ax4.set_ylim([X.min(), ca])
#ax1.set_aspect("equal")
cbar = fig.colorbar(cmesh, ax=ax4)
cbar.set_label("Velocity (m/s)")
ax4.set_xlabel("Parallel to B (m)")
ax4.set_ylabel("Radial (m)")

# Plot of the C density 1 cm above the limiter.
#ax5.plot(y_coords[:-1], above_ddlim)
#ax5.set_xlabel("Parallel to B (m)")
#ax5.set_ylabel("C Density (m-3)")
#ax5.set_xlim([-cl, cl])
#ax5.set_xlim([-3, 3])
#ax5.set_yscale("log")

cmesh = ax5.pcolormesh(Y, X, Z, cmap="inferno", norm=norm)
ax5.fill_between(-qedges1, qxs, X.min(), color="grey")
ax5.fill_between(qedges2, qxs, X.min(), color="grey")
ax5.plot(-qedges1, qxs, color="grey", lw=2)
ax5.plot(qedges2, qxs, color="grey", lw=2)
ax5.set_xlim([-10, 10])
ax5.set_ylim([X.min(), ca])
#if equal:
#    ax1.set_aspect("equal")
cbar = fig.colorbar(cmesh, ax=ax5)
cbar.set_label("C Density (m-3)")
ax5.set_xlabel("Parallel to B (m)")
ax5.set_ylabel("Radial (m)")

# Plot of erosion along the limiter surface. NERODS3 may be in units of particles/m2.
#cmesh = ax6.pcolormesh(O, P, nerods3_masked*absfac, cmap="inferno")
#cbar = fig.colorbar(cmesh, ax=ax6)
#cbar.set_label("Erosion particles/m2/s?")
#ax6.set_aspect("equal")

# Plot of impact energy and the flux.

fig.tight_layout()
fig.show()
