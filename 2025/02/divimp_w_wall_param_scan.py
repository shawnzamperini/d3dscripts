import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# Figure of merit. This is what is used to compare the different scenarios
# to each other. Options are: avg_conc
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
dperps = [0.1, 0.3, 0.6, 1.0, 5.0]
ne_charges = [3, 5, 7, 9, 11]
ne_fluxes = [0.005, 0.01, 0.03, 0.05, 0.10]

loop_count = 1
for dperp in dperps:
	if loop_count > 80: break
	for ne_charge in ne_charges:
		if loop_count > 80: break
		for ne_flux in ne_fluxes:
			if loop_count > 80: break

			# Load case
			ncpath = "/home/zamp/oedge_files/d3d-w-wall-param-scan-{}.nc".format(loop_count)
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
				ax.plot(s, cw_s, lw=3, alpha=0.3)

				# We may want to ignore the ends of the flux tubes
				s_keep = np.logical_and(s > s_trim, s < (s.max() - s_trim))
				fom_vals.append(cw_s[s_keep].mean())
				fom_stds.append(cw_s[s_keep].std())
			
			# Figure of merit: Average W density
			elif fom == "avg_nw":
				ax.plot(s, nw_s, lw=3, alpha=0.3)

				# We may want to ignore the ends of the flux tubes
				s_keep = np.logical_and(s > s_trim, s < (s.max() - s_trim))
				fom_vals.append(nw_s[s_keep].mean())
				fom_stds.append(nw_s[s_keep].std())

			# Record the values of each parameter
			xvalues.append([dperp, ne_charge, ne_flux])

			# Increment
			loop_count += 1

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

ax.legend()
fig.tight_layout()
fig.show()


# Now that we have a large enough set of values, let's perform a multivariable
# curve fit to derive a scaling law
def model(x, a, b, c):
	
	# These are: D, T, Z and c
	x1, x2, x3 = x

	# Concentration is by definition 1.0 exponent here, so effectively remove it
	# as a free parameter.
	return x1**a * x2**b * x3**c

popt, pcov = curve_fit(model, xvalues.T, fom_vals,
	bounds=((-10, -10, -10), (10, 10, 10)))
fom_fit = model(xvalues.T, *popt)

print("Fit Parameters")
print("  dperp:     {:.2f}".format(popt[0]))
print("  ne_charge: {:.2f}".format(popt[1]))
print("  ne_flux:   {:.2f}".format(popt[2]))

# Still growing this plot, but it's the figure of merit plot
fig, ax = plt.subplots()
ax.scatter(fom_fit, fom_vals, s=75, color="tab:red", edgecolors="k")
ax.plot([fom_fit.min(), fom_fit.max()], [fom_vals.min(), fom_vals.max()], 
	lw=3, color="k", linestyle="--") 
ax.set_xlabel("Fit", fontsize=fontsize)
ax.grid(alpha=0.3)
ax.set_ylabel(fom, fontsize=fontsize)
ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_xlim([1e-4, 1e-2])
#ax.set_ylim([1e-4, 1e-2])
fig.tight_layout()
fig.show()

