# Simple script to combine some repeat DIVIMP runs and plot them.
import oedge_plots
import netCDF4


# Load each as an OedgePlots object.
root = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-pdf10b-ring-"
#root = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-corr-"
ops = []
for letter in ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]:
    ops.append(oedge_plots.OedgePlots("{}{}.nc".format(root, letter)))

# Pull out the total W densities and combine into one.
tot_dens = ops[0].read_data_2d("DDLIMS", charge="all")
for i in range(1, len(ops)):
    tot_dens += ops[i].read_data_2d("DDLIMS", charge="all")

# Average values.
tot_dens = tot_dens / len(ops)

# ABSFAC is this.
absfac = 7.89e17
tot_dens = tot_dens * absfac
vmin = 1e-7 * absfac
vmax = 1e-3 * absfac

# Plotting. Total particles = 6x50,000 + 6x100,000 = 900,000 (too many!). 
ops[0].plot_contour_polygon("KTEBS", own_data=tot_dens, vmin=vmin, vmax=vmax,
    cmap="nipy_spectral", cbar_label="W Density (m-3)", normtype="log")
