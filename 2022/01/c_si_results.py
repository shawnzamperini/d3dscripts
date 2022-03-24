# This script is to be used to combine two indentical 3DLIM runs differing just
# in the one uses physical sputtering and the other uses chemical.
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm, SymLogNorm, BoundaryNorm
from matplotlib.lines import Line2D


title = "Best Case Scenario: RCP input as-is, Ti = Te\n3% C background"
#title = "Worst Case Scenario: RCP input shifted 2 cm, Ti = 3Te\n3% C background"

# Load nc objects here.
# SiC - Best case scenario: 006, 008, 007 and 012, 013
# SiC - Worst case scenario: 009, 011, 010 and 014, 015
# Graphite - Best case scenario: 012, 013, dummy and dummy, dummy
# Graphite - Worst case scenario: 014, 015, dummy and dummy, dummy
ncpath_phys = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-d-sic-c-006.nc"
ncpath_chem = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-d-sic-c-008.nc"
ncpath_sili = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-d-sic-c-007.nc"
ncpath_grap_phys = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-d-sic-c-012.nc"
ncpath_grap_chem = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-d-sic-c-013.nc"
nc_phys = netCDF4.Dataset(ncpath_phys)
nc_chem = netCDF4.Dataset(ncpath_chem)
nc_sili = netCDF4.Dataset(ncpath_sili)
nc_grap_phys = netCDF4.Dataset(ncpath_grap_phys)
nc_grap_chem = netCDF4.Dataset(ncpath_grap_chem)

def load_output_dict(nc, dist_above_lim=0.01, include_neutrals=True,
  pol_idx=21):
    """
    Load a dictionary with a select portion of the netcdf data, as well as
    some data that has been post processed to pull out trends.

    nc (Dataset): The netCDF file.
    dist_above_lim (float): Distance above limiter to provide parallel DDLIM3
      data for.
    include_neutrals (bool): Include neutrals in DDLIM3 data.
    pol_idx (int): Poloidal index of which to pull data from.
    """

    output = {}

    # Extract some data up front.
    xs = nc.variables["XS"][:].data
    ys = nc.variables["YS"][:].data
    xwids = nc.variables["XWIDS"][:].data
    ywids = nc.variables["YWIDS"][:].data
    youts = nc.variables['YOUTS'][:].data
    xouts = nc.variables['XOUTS'][:].data
    absfac = nc.variables["ABSFAC"][:].data
    qtembs = nc.variables["QTEMBS"][:][0][:-2]  # Last 2 values seem to be junk.
    qtembsi = nc.variables["QTEMBSI"][:][0][:-2]
    cthetb = float(nc.variables["CTHETB"][:].data)
    qrnbs = nc.variables["QRNBS"][:][0][:-2]
    odouts = nc.variables["ODOUTS"][:]
    cl = float(nc['CL'][:].data)
    ca  = float(nc['CA'][:].data)
    caw = float(nc['CAW'][:].data)

    # The limiter locations.
    qxs = nc.variables["QXS"][:].data
    qedges1 = nc.variables["QEDGES"][:].data[0]
    qedges2 = nc.variables["QEDGES"][:].data[1]

    # Erosion data. 0 = deposition. 1 = primary erosion.
    # 2 = erosion + self-sputtering. 4 = net erosion
    nerods3 = nc.variables["NERODS3"][:].data[2][pol_idx] * absfac

    # Only need the values up to len(qedges). Trim zeros. Want one zero at the
    # end for fill between.
    nonzero = np.nonzero(qedges1)
    qxs = qxs[:len(qedges1)][nonzero]
    qedges1 = qedges1[nonzero]
    qedges2 = qedges2[nonzero]
    qxs = np.append(qxs, 0)
    qedges1 = np.append(qedges1, 0)
    qedges2 = np.append(qedges2, 0)

    # Mirror the Y data around 0.
    ys = np.append(np.append(-ys[::-1], 0), ys)
    ywids = np.append(np.append(-ywids[::-1], 0), ywids)
    youts = np.append(np.append(-youts[::-1], 0), youts)

    # The longest part of the process is summing ddlim3 across the charge states.
    if include_neutrals:
        ddlim3 = nc.variables["DDLIM3"][:].data[pol_idx, 1:, :, :].sum(axis=0)
    else:
        ddlim3 = nc.variables["DDLIM3"][:].data[pol_idx, 2:, :, :].sum(axis=0)
    ddlim3 = np.ma.masked_where(ddlim3<=0, ddlim3)

    # Scale to m-3.
    ddlim3 = ddlim3 * absfac

    # We then need to trim all the trailing zeros on the arrays.
    xkeep_min = np.nonzero(xs)[0].min()
    xkeep_max = np.nonzero(xs)[0].max()
    ykeep_min = np.nonzero(ys)[0].min()
    ykeep_max = np.nonzero(ys)[0].max()
    ys     = ys[ykeep_min:ykeep_max]
    ywids  = ywids[ykeep_min:ykeep_max]
    youts  = youts[ykeep_min:ykeep_max]
    ddlim3 = ddlim3[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

    # The net erosion at the poloidal index. The units are atoms/m2/s, where
    # the area is the surface area (not parallel or anything).
    nerods = nerods3[pol_idx] * absfac

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
    ddlim3 = ddlim3[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

    # X, Y 2D arrays fro the plots.
    x_coords = np.append(xouts, xouts[-1] + xwids[-1] / 2)
    y_coords = np.append(youts, youts[-1] + ywids[-1] / 2)
    X, Y = np.meshgrid(x_coords, y_coords)

    # Create Y coordinates for the DDLIM3 data above the limiter.
    above_idx = np.argmin(np.abs(xouts - dist_above_lim))
    above_ddlim3 = ddlim3[:, above_idx]

    # Calculate the surface flux as calculated in the code. Need to mirror it
    # so it aligns with odouts and makes sense.
    cs = 9.79e3 * np.sqrt((qtembs + qtembsi) / 2)
    csintb = np.sin(np.radians(cthetb))
    flux = qrnbs * cs * csintb
    half = int(len(flux)/2)
    flux = np.append(flux[half:], flux[half:][::-1])

    output["X"] = X
    output["Y"] = Y
    output["ddlim3"] = ddlim3
    output["above_y"] = y_coords
    output["above_ddlim3"] = above_ddlim3
    output["odouts"] = odouts
    output["nerods3"] = nerods3
    output["qtembs"] = qtembs
    output["qtembsi"] = qtembsi
    output["qrnbs"] = qrnbs
    output["qedges1"] = qedges1
    output["qedges2"] = qedges2
    output["qxs"] = qxs
    output["ca"] = ca
    output["caw"] = caw
    output["cl"] = cl
    output["flux"] = flux

    return output

print("Loading physical results...")
phys = load_output_dict(nc_phys)
print("Loading chemical results...")
chem = load_output_dict(nc_chem)
print("Loading silicon results...")
sili = load_output_dict(nc_sili)
outputs = [phys, chem, sili]

# Create figure with large axes behind everything for common labels.
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.set_xlabel("Parallel (m)", fontsize=14)
ax.set_ylabel("Radial (m)\n", fontsize=14)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

# The actual plots.
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8,5))

# Setup for plotting loop.
norm = LogNorm(vmin=1e14, vmax=1e18)
axs = [ax1, ax2, ax3, ax4]
xlim1 = [-0.15, 0.15]; equal = True
#xlim1 = [-45, 45]; equal = False
ncolors = 10
#bounds = np.power(10.0, np.arange(14, 19))
bounds = np.geomspace(1e14, 1e18, ncolors)
cmap = plt.cm.get_cmap('inferno', ncolors)
norm = BoundaryNorm(boundaries=bounds, ncolors=ncolors)
yticks = np.power(10.0, np.arange(14, 19))
yticklabels = [f'$10^{{{np.log10(b):.0f}}}$' for b in yticks]

for i in range(0, 4):

    # Sum of C physical + chemical.
    if i == 3:
        X = outputs[0]["X"]
        Y = outputs[0]["Y"]
        Z = outputs[0]["ddlim3"].data + outputs[1]["ddlim3"].data
        Z = np.ma.masked_where(Z<=0, Z)
        qedges1 = outputs[0]["qedges1"]
        qedges2 = outputs[0]["qedges2"]
        qxs = outputs[0]["qxs"]
        ca = outputs[0]["ca"]
        Z_sic = Z
    else:
        X = outputs[i]["X"]
        Y = outputs[i]["Y"]
        Z = outputs[i]["ddlim3"]
        qedges1 = outputs[i]["qedges1"]
        qedges2 = outputs[i]["qedges2"]
        qxs = outputs[i]["qxs"]
        ca = outputs[i]["ca"]

    cmesh = axs[i].pcolormesh(Y, X, Z, cmap=cmap, norm=norm)
    axs[i].fill_between(-qedges1, qxs, X.min(), color="grey")
    axs[i].fill_between(qedges2, qxs, X.min(), color="grey")
    axs[i].plot(-qedges1, qxs, color="grey", lw=2)
    axs[i].plot(qedges2, qxs, color="grey", lw=2)
    axs[i].set_xlim(xlim1)
    axs[i].set_ylim([X.min(), ca])
    if equal:
        axs[i].set_aspect("equal")

# Some fine tuning around the subplots.
ax2.set_yticks([])
ax4.set_yticks([])
ax1.set_title("C Physical Sputtering")
ax2.set_title("C Chemical Sputtering")
ax3.set_title("Si Physical Sputtering")
ax4.set_title("C Physical + Chemical")

# https://stackoverflow.com/questions/13784201/how-to-have-one-colorbar-for-all-subplots
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(cmesh, cax=cbar_ax)
cbar_ax.set_ylabel("Density (m-3)", fontsize=14)
cbar_ax.set_yticks(yticks)
cbar_ax.set_yticklabels(yticklabels)
fig.suptitle(title)
#fig.tight_layout()
fig.show()

# Bonus plot: Comparison to pure graphite limiter case.
print("Loading graphite physical...")
grap_phys = load_output_dict(nc_grap_phys)
print("Loading graphite chemical...")
grap_chem = load_output_dict(nc_grap_chem)

Z_grap = grap_phys["ddlim3"].data + grap_chem["ddlim3"].data
Z_grap = np.ma.masked_where(Z_grap<=0, Z_grap)
Z_ratio = ((Z_sic - Z_grap) / Z_grap) * 100

ncolors = 11
bounds = np.linspace(-100, 100, ncolors)
cmap = plt.cm.get_cmap('coolwarm', ncolors)
norm = BoundaryNorm(boundaries=bounds, ncolors=ncolors)

fig, ax = plt.subplots()
cmesh = ax.pcolormesh(Y, X, Z_ratio, cmap=cmap, norm=norm)
ax.fill_between(-qedges1, qxs, X.min(), color="grey")
ax.fill_between(qedges2, qxs, X.min(), color="grey")
ax.plot(-qedges1, qxs, color="grey", lw=2)
ax.plot(qedges2, qxs, color="grey", lw=2)
ax.set_xlim(xlim1)
ax.set_ylim([X.min(), ca])
if equal:
    ax.set_aspect("equal")
cbar = fig.colorbar(cmesh, ax=ax, location="top")
cbar.set_label("Change in C Density (%)")
ax.set_xlabel("Parallel to B (m)")
ax.set_ylabel("Radial (m)")
fig.tight_layout()
fig.show()

# Another bonus plot. Comparison of the erosion profiles for SiC to graphite.
odouts = outputs[0]["odouts"]
#nerods_tot_sic = outputs[0]["nerods3"] + outputs[1]["nerods3"]
#nerods_tot_gra = grap_phys["nerods3"] + grap_chem["nerods3"]
nerods_tot_sic = outputs[0]["nerods3"]
nerods_tot_gra = grap_phys["nerods3"]
nerods_sic_ratio = nerods_tot_sic / outputs[0]["flux"]
nerods_gra_ratio = nerods_tot_gra / grap_phys["flux"]

fig, ax = plt.subplots()

ax.plot(odouts, nerods_sic_ratio, label="SiC", color="tab:purple", lw=2)
ax.plot(odouts, nerods_gra_ratio, label="Graphite", color="tab:red", lw=2)
ax.axhline(0, color="k", lw=2, linestyle="--")
ax.legend(fontsize=12)
ax.set_xlabel("Distance along limiter from tip (m)", fontsize=14)
ax.set_ylabel("Gross Physical Erosion / D Flux", fontsize=14)
ax.grid()
#ax.set_ylim(-1e22, 0.3e22)
ax.set_ylim([0, None])

fig.tight_layout()
fig.show()
