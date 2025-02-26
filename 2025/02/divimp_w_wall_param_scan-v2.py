import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# Figure of merit. This is what is used to compare the different scenarios
# to each other. Options are: avg_conc, avg_nw
fom = "avg_conc"

# Load base case
#ncpath = "/home/zamp/oedge_files/d3d-w-wall-dperp-scan-001.nc"
#op = oedge_plots.OedgePlots(ncpath)

# Load data so we can manipulate it before plotting
#nw = op.read_data_2d("DDLIMS", charge="all", scaling=op.absfac)
#rings = op.read_data_2d("KTEBS", scaling="Ring")

# We want to plot the W density, but mask the data where the flows are 
# that prevent transport from really occuring (rings less than 40, non-
# inclusive).
ring = 18
#nw[np.logical_or(rings < ring, rings > 158)] = np.nan
#op.plot_contour_polygon("KTEBS", own_data=nw, normtype="log", 
#	cbar_label="W Density (m-3)", vmin=1e13, vmax=1e16)

# Create figure for the plotting below
fontsize = 16
fig, ax = plt.subplots()

# Lists to hold our figures of merit as we extract them from each case
fom_vals = []
fom_stds = []
xvalues = []

# The parameters that were scanned
dperps = [3.907, 1.081, 1.510, 2.956, 0.305, 1.241, 1.782, 0.550, 3.039, 0.348, 4.443, 4.897, 0.711, 4.543, 4.447, 2.874, 0.616, 4.368, 3.565, 2.515, 3.946, 0.248, 1.255, 4.218, 3.266, 0.171, 2.568, 3.900, 4.616, 3.087, 4.554, 3.258, 1.972, 3.642, 4.860, 2.988, 4.867, 0.864, 2.244, 2.061, 1.179, 4.612, 3.009, 2.094, 4.173, 1.155, 4.286, 2.353, 0.177, 1.634, 2.149, 4.605, 4.770, 3.941, 0.484, 1.368, 2.661, 1.653, 3.202, 1.459, 4.471, 4.476, 3.691, 0.840, 3.215, 4.186, 3.834, 2.290, 4.071, 2.930, 1.684, 4.792, 2.401, 4.609, 4.981, 1.669, 3.130, 0.346, 2.826, 1.575, 1.120, 4.577, 4.396, 3.034, 1.777, 4.540, 1.069, 2.206, 0.796, 1.737, 4.612, 4.823, 0.246, 1.113, 2.755, 4.220, 0.337, 4.612, 0.301, 3.772]
ne_charges = [7, 2, 7, 3, 5, 3, 3, 6, 3, 6, 2, 8, 5, 9, 6, 8, 2, 4, 3, 6, 5, 2, 5, 7, 6, 7, 6, 5, 8, 6, 2, 2, 9, 4, 7, 6, 9, 8, 8, 4, 4, 8, 8, 4, 9, 9, 6, 8, 5, 6, 7, 6, 4, 3, 7, 2, 2, 4, 5, 4, 2, 5, 3, 4, 3, 7, 5, 4, 7, 9, 7, 7, 3, 6, 8, 9, 3, 8, 6, 6, 5, 2, 4, 9, 6, 6, 9, 5, 6, 7, 6, 5, 2, 5, 4, 5, 7, 6, 7, 5]
ne_fluxes = [0.082, 0.037, 0.019, 0.064, 0.021, 0.044, 0.073, 0.014, 0.099, 0.023, 0.009, 0.012, 0.048, 0.071, 0.027, 0.079, 0.058, 0.056, 0.040, 0.012, 0.022, 0.053, 0.030, 0.028, 0.070, 0.029, 0.064, 0.051, 0.029, 0.058, 0.066, 0.093, 0.085, 0.008, 0.040, 0.071, 0.077, 0.038, 0.014, 0.069, 0.071, 0.075, 0.064, 0.059, 0.064, 0.042, 0.018, 0.091, 0.051, 0.056, 0.021, 0.037, 0.081, 0.014, 0.059, 0.011, 0.029, 0.075, 0.076, 0.051, 0.076, 0.007, 0.095, 0.024, 0.027, 0.049, 0.022, 0.011, 0.073, 0.059, 0.077, 0.062, 0.058, 0.020, 0.058, 0.042, 0.066, 0.087, 0.010, 0.091, 0.022, 0.022, 0.013, 0.070, 0.025, 0.016, 0.026, 0.095, 0.067, 0.018, 0.087, 0.064, 0.061, 0.095, 0.043, 0.021, 0.097, 0.064, 0.074, 0.027]

for i in range(0, len(dperps)):

	# Load case
	#ncpath = "/home/zamp/oedge_files/d3d-w-wall-param-scan-v2-{}.nc".format(i+1)
	ncpath = "/home/zamp/oedge_files/d3d-w-wall-param-scan-v3-all-{}.nc".format(i+1)
	print(ncpath)
	op = oedge_plots.OedgePlots(ncpath)

	# Load needed density to derive figure of merit at desired ring
	s, nw_s = op.along_ring(ring, "DDLIMS", charge="all", 
		scaling=op.absfac, plot_it=False)
	s, ne_s = op.along_ring(ring, "KNBS", plot_it=False)
	cw_s = nw_s / (nw_s + ne_s)

	# We may want to ignore the ends of the flux tubes
	s_trim = 0
	
	# Figure of merit: Average W concentration
	if fom == "avg_conc":
		ax.plot(s, cw_s, lw=4, alpha=0.3, color="k")
		ax.plot(s, cw_s, lw=3, alpha=0.3, color="tab:red")

		# We may want to ignore the ends of the flux tubes
		s_keep = np.logical_and(s > s_trim, s < (s.max() - s_trim))
		fom_vals.append(cw_s[s_keep].mean())
		fom_stds.append(cw_s[s_keep].std())
	
	# Figure of merit: Average W density
	elif fom == "avg_nw":
		ax.plot(s, nw_s, lw=3, alpha=0.3, color="tab:red")

		# We may want to ignore the ends of the flux tubes
		s_keep = np.logical_and(s > s_trim, s < (s.max() - s_trim))
		fom_vals.append(nw_s[s_keep].mean())
		fom_stds.append(nw_s[s_keep].std())

	# Record the values of each parameter
	xvalues.append([dperps[i], ne_charges[i], ne_fluxes[i]])

fom_vals = np.array(fom_vals)
fom_stds = np.array(fom_stds)
xvalues = np.array(xvalues)

#ax.set_yscale("log")
ax.set_xlabel("Distance from inner target (m)", fontsize=fontsize)

# Format y axis based off figure of merit.
if fom == "avg_conc":
	ax.set_ylabel("W Concentration", fontsize=fontsize)
	#ax.set_ylim([0, 1e-5])
elif fom == "avg_nw":
	ax.set_ylabel("W Density (m-3)", fontsize=fontsize)

#ax.legend()
ax.set_yscale("log")
#ax.set_ylim([5e-6, 3e-2])
ax.set_ylim([1e-4, 1e0])
fig.tight_layout()
fig.show()


# Now that we have a large enough set of values, let's perform a multivariable
# curve fit to derive a scaling law
def model(x, a, b, c, d):
	
	# These are: D, T, Z and c
	x1, x2, x3 = x

	# Concentration is by definition 1.0 exponent here, so effectively remove it
	# as a free parameter.
	return x1**a * x2**b * x3**c * d

popt, pcov = curve_fit(model, xvalues.T, fom_vals,
	bounds=((-10, -10, -10, -np.inf), (10, 10, 10, np.inf)))
fom_fit = model(xvalues.T, *popt)

print("Fit Parameters")
print("  dperp:     {:.2f}".format(popt[0]))
print("  ne_charge: {:.2f}".format(popt[1]))
print("  ne_flux:   {:.2f}".format(popt[2]))

# Create a text object to place on the plot of the fit
if fom == "avg_conc":
	lhs = r"$\mathdefault{c_W}$"
elif fom == "avg_nw":
	lhs = r"$\mathdefault{n_W}$"
text = lhs + r"$\mathdefault{\sim D_r^{" + "{:.1f}".format(popt[0]) + r"}\ Z_{Ne}^{" + "{:.1f}".format(popt[1]) + r"}\ f_{Ne}^{" + "{:.1f}".format(popt[2]) + r"}}$"

# Still growing this plot, but it's the figure of merit plot
fontsize = 18
fig, ax = plt.subplots()
ax.scatter(fom_fit, fom_vals, s=75, color="tab:red", edgecolors="k")
ax.plot([1e-4, 3e-1], [1e-4, 3e-1], 
	lw=3, color="k", linestyle="--")
ax.set_xlabel(text, fontsize=fontsize)
ax.set_ylabel(lhs, fontsize=fontsize)
ax.tick_params(axis="both", labelsize=fontsize-4)
ax.grid(alpha=0.3)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim([1e-4, 3e-1])
ax.set_ylim([1e-4, 3e-1])
fig.tight_layout()
fig.show()

