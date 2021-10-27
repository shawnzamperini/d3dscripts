# Copy/paste with modifications of 2021/09/ddlim3_gif.py.
# This script is designed to create a gif of progressing through the DDLIM3
# array from the "top-down". This is to show how the parallel-perpendicular
# impurity distribution changes as one progresses radially outward.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches
import netCDF4
import imageio
import os


# Provide inputs and plot options here.
ncpath = "/Users/zamperini/Documents/lim_runs/actual-184527-004.nc"  # This is a simulation using connection lengths at toroidal angle = 0.
#ncpath = "/Users/zamperini/Documents/lim_runs/two_ramps_test.nc"
cmap = "magma"
rad_sum_width = 0.005  # If None just go one radial bin at a time.
#rad_sum_width = None
xlims = [-51, 51]  # The parallel limits.
xlabels = np.arange(-51, 51, 20, dtype=int)  # Labels for parallel axis.
ylims = [2.26, 2.40]  # The radial limits.
extra_frames = 1  # Additional frames to interpolate between frames.
cp_width = 0.5  # Width of the "CP" surface.
r_tip = 2.29   # The R value of the origin to convert to easier coordinates.
p_to_z = [0.97636342, -0.188]  # Z = [0] * pbin + [1], i.e. linear fit to go
                               # from P to Z as figured out in
                               # bounds_file_from_mafot.py.
plot_vel = True

# Load the NetCDF and pull out various arrays.
print("Loading netCDF file and extracting variables...")
nc = netCDF4.Dataset(ncpath)
ps = nc.variables["PS"][:].data
xs = nc.variables["XS"][:].data
ys = nc.variables["YS"][:].data
ywids = nc.variables["YWIDS"][:].data
xwids = nc.variables["XWIDS"][:].data
pwids = nc.variables["PWIDS"][:].data
ddlim3 = nc.variables["DDLIM3"][:].data
vp = nc.variables["velplasma_4d_1"][:].data

# Warning if the width of the sum range is smaller than the smallest radial
# bin width...
min_xwid = xwids[np.nonzero(xwids)].min()

# Sum across all charge states.
print("Summing across charge states...")
ddlim3 = ddlim3.sum(axis=1)

# Calculate centers of bins. Special treatment for the Y coordinate since it
# needs to be mirrored and a zero added.
rad_locs = r_tip - (xs - xwids / 2)
pol_locs = p_to_z[0] * (ps - pwids / 2) + p_to_z[1]  # This is actually Z now.
tmp = ys - ywids / 2
par_locs = np.append(np.append(-tmp[::-1], 0), tmp)
all_par_locs = par_locs
all_rad_locs = (xs - xwids / 2)
all_pol_locs = (ps - pwids / 2)

# Cut down on the array sizes by dropping unecessary zeros. The +1 is so that we
# include the max value in the indexing.
y_keep_start = np.nonzero(par_locs)[0].min()
y_keep_end   = np.nonzero(par_locs)[0].max() + 1
x_keep_start = np.nonzero(rad_locs)[0].min()
x_keep_end   = np.nonzero(rad_locs)[0].max() + 1
p_keep_start = np.nonzero(pol_locs)[0].min()
p_keep_end   = np.nonzero(pol_locs)[0].max() + 1
rad_locs = rad_locs[x_keep_start:x_keep_end]
pol_locs = pol_locs[p_keep_start:p_keep_end]
par_locs = par_locs[y_keep_start:y_keep_end]
ddlim3 = ddlim3[p_keep_start:p_keep_end, y_keep_start:y_keep_end,
  x_keep_start:x_keep_end]
vp = vp[y_keep_start:y_keep_end, x_keep_start:x_keep_end,
  p_keep_start:p_keep_end]

# A plot in the R, Z plane of the impurity density. Only really makes sense
# I think to do this at the origin here since that's the only R, Z plane that
# is I guess in the normal 2D shape.

# rad_locs contains a lot of trailing zeros (or actually r_origin since we do
# r_origin - rbin), and we need to ditch those.
start_ignore = 49
tmp_rad_locs = rad_locs[:start_ignore]

origin_idx = np.where(par_locs==0)[0][0] + 1
ddlim3_origin = ddlim3[:, origin_idx, :start_ignore]
vp_ex = vp[:, :start_ignore, 10]

fig, ax = plt.subplots()

# pol_locs actually Z if you see above, just lazy.
#cont = ax.pcolormesh(tmp_rad_locs, pol_locs, ddlim3_origin, shading="auto")
cont = ax.pcolormesh(tmp_rad_locs, par_locs, vp_ex, shading="auto", vmin=-15000,
  vmax=15000, cmap="coolwarm")
cbar = fig.colorbar(cont)
fig.tight_layout()
fig.show()


# Create bins of each radial range we will sum ddlim3 over, and assign which
# indices in the rdial dimension belong to which bin.
summed_ddlim3 = []; avg_r = []
if rad_sum_width != None:
    xbins_idx = []
    xbins = np.arange(rad_locs.min(), rad_locs.max(), rad_sum_width)
    for i in range(0, len(xbins)):
        if i == len(xbins)-1:
            idx = np.where(rad_locs >= xbins[i])
            avg_r.append(xbins[i])
        else:
            idx = np.where(np.logical_and(rad_locs >= xbins[i],
              rad_locs < xbins[i+1]))
            avg_r.append((xbins[i] + xbins[i+1]) / 2)
        xbins_idx.append(idx[0])

    # Performing the summing for each radial bin.
    for idx in xbins_idx:
        summed_ddlim3.append(ddlim3[:,:,idx].sum(axis=2))

else:
    xbins_idx = np.arange(0, ddlim3.shape[2], dtype=int)
    for idx in xbins_idx:

        # No sum if we just want this radial bin's values.
        summed_ddlim3.append(ddlim3[:,:,idx])

# Reverse so it goes from top to bottom (order of X bins go bottom to top) and
# convert to numpy array.
#summed_ddlim3.reverse()
#avg_r.reverse()
summed_ddlim3 = np.array(summed_ddlim3)

# Grab some values for the plot properties.
#vmin = summed_ddlim3[np.nonzero(summed_ddlim3)].min()
vmin = 0.00001
vmax = summed_ddlim3.max()
#lev_exp = np.linspace(np.floor(np.log10(vmin)-1),
#  np.ceil(np.log10(vmax)+1), 10)
#levs = np.power(10, lev_exp)
levs = np.geomspace(vmin, vmax, 10)
X1, Y1 = np.meshgrid(par_locs, pol_locs)
X2, Y2 = np.meshgrid(par_locs, rad_locs)

# Z2 is the summed density across the poloidal range. More or less used here
# qualitatively.
Z2 = ddlim3.sum(axis=0)
Z2 = np.ma.masked_where(Z2 <= 0, Z2).T
levs2 = np.geomspace(0.0001, Z2.max(), 15)

# Now create a list of plots and save them for putting into a gif later.
print("Generating plots...")
fnames = []
for i in range(0, len(summed_ddlim3)-1):
    print(" plot {:}...".format(i))

    # This loops calculate the difference in ddlim3 for each frame location
    # and then interpolate extra_frame between them to slow/smooth the gif out.
    ddlim3_diff = summed_ddlim3[i+1] - summed_ddlim3[i]
    r_diff = avg_r[i+1] - avg_r[i]
    for j in range(0, extra_frames+1):

        # Calculate the interpolated frame and plot it. Mask to prevent a
        # warning with log stuff.
        interp_ddlim3 = summed_ddlim3[i] + (ddlim3_diff / extra_frames) * j
        interp_ddlim3 = np.ma.masked_where(interp_ddlim3 <= 0, interp_ddlim3)
        interp_r = avg_r[i] + (r_diff / extra_frames) * j

        # Plot the 2D distribution.
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        c = ax1.contourf(X1, Y1, interp_ddlim3, levs, cmap=cmap,
          norm=colors.LogNorm())
        #c = ax1.pcolormesh(X1, Y1, interp_ddlim3, cmap=cmap,
        #  norm=colors.LogNorm(vmin=vmin, vmax=vmax), shading="auto")

        # Add in a little CP as well.
        cp1 = patches.Rectangle((-0.05, -cp_width), 0.1, cp_width*2,
          color="grey", ec="k")
        #ax1.add_patch(cp1)

        # Miscellanous plot details.
        ax1.set_facecolor("grey")
        #cbar = fig.colorbar(c, ax=ax1)
        #cbar.set_label("C Density (arbitrary)")
        #cbar.set_ticks([])  # Arbitrary units anyways.
        ax1.set_xlim(xlims)
        ax1.set_ylim([-0.7, 0.3])
        ax1.set_xlabel("Distance from floor (m)")
        ax1.set_ylabel("Z (m)")
        ax1.set_xticks(xlabels)
        ax1.set_xticklabels(xlabels)
        ax1.set_title("Top-down view")

        # This is an accompanying side plot to show the radial position of each
        # slice.
        ax2.set_facecolor("grey")
        ax2.contourf(X2, Y2, Z2, levs2, norm=colors.LogNorm(), cmap=cmap)
        ax2.set_xlabel("Distance from floor (m)")
        ax2.set_ylabel("R (m)")
        cp2 = patches.Rectangle((-0.05, -0.15), 0.1, 0.15, color="grey", ec="k")
        #ax2.add_patch(cp2)
        ax2.spines["top"].set_visible(False)
        ax2.spines["bottom"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        #ax2.tick_params(axis="both", which="both", bottom=False, top=False,
        #  left=False, right=False, labelbottom=False, labelleft=False)
        ax2.set_xlim(xlims)
        ax2.set_ylim(ylims)
        ax2.axhline(interp_r, color="k", lw=3)
        ax2.axhline(interp_r, color="r", lw=2)
        ax2.set_xticks(xlabels)
        ax2.set_xticklabels(xlabels)
        ax2.set_title("Side view")
        ax2.invert_yaxis()

        # Save figure and close.
        fname = "ddlim3_plots/{:}_{:}.png".format(i,j)
        fnames.append(fname)
        fig.savefig(fname)
        plt.close(fig)

# Now build the gif and clean up by removing the plots.
print("Saving as gif...")
with imageio.get_writer('ddlim3.gif', mode="I") as writer:
    for fname in fnames:
        image = imageio.imread(fname)
        writer.append_data(image)
for fname in fnames:
    os.remove(fname)
