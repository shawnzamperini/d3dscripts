import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# Figure of merit. This is what is used to compare the different scenarios
# to each other. Options are: avg_conc
fom = "avg_conc"

# Load base case
ncpath = "/home/zamp/oedge_files/d3d-w-wall-dperp-scan-001.nc"
op = oedge_plots.OedgePlots(ncpath)

# Load data so we can manipulate it before plotting
nw = op.read_data_2d("DDLIMS", charge="all", scaling=op.absfac)
rings = op.read_data_2d("KTEBS", scaling="Ring")

# We want to plot the W density, but mask the data where the flows are 
# that prevent transport from really occuring (rings less than 40, non-
# inclusive).
ring = 41
nw[np.logical_or(rings < ring, rings > 158)] = np.nan
op.plot_contour_polygon("KTEBS", own_data=nw, normtype="log", 
	cbar_label="W Density (m-3)", vmin=1e13, vmax=1e16)

# Create figure for the plotting below
fontsize = 16
fig, ax = plt.subplots()

fom_vals = []
fom_stds = []
dperp = [0.3, 0.6, 1.0, 5.0, 0.1]
for i in range(1, 6):
	ncpath = "/home/zamp/oedge_files/d3d-w-wall-dperp-scan-00{}.nc".format(i)
	op = oedge_plots.OedgePlots(ncpath)

	# Then for comparisons, we want to look at the W conc at the furthest in
	# ring we care about (40)
	s, nw_s = op.along_ring(ring, "DDLIMS", charge="all", scaling=op.absfac, 
		plot_it=False)
	s, ne_s = op.along_ring(ring, "KNBS", plot_it=False)
	cw_s = nw_s / (nw_s + ne_s)

	if fom == "avg_conc":
		ax.plot(s, cw_s, lw=3, color="k")
		ax.plot(s, cw_s, label="D = {} m2/s".format(dperp[i-1]), lw=2)

		# For average concentration, let's ignore the first and last 10m of the
		# flux tube
		s_trim = 10
		s_keep = np.logical_and(s > s_trim, s < (s.max() - s_trim))
		fom_vals.append(cw_s[s_keep].mean())
		fom_stds.append(cw_s[s_keep].std())
	
	elif fom == "avg_nw":
		ax.plot(s, nw_s, lw=3, color="k")
		ax.plot(s, nw_s, label="D = {} m2/s".format(dperp[i-1]), lw=2)

		# For average concentration, let's ignore the first and last 10m of the
		# flux tube
		s_trim = 10
		s_keep = np.logical_and(s > s_trim, s < (s.max() - s_trim))
		fom_vals.append(nw_s[s_keep].mean())
		fom_stds.append(nw_s[s_keep].std())
fom_vals = np.array(fom_vals)
fom_stds = np.array(fom_stds)

#ax.set_yscale("log")
ax.set_xlabel("Distance from inner target (m)", fontsize=fontsize)

# Format y axis based off figure of merit.
if fom == "avg_conc":
	ax.set_ylabel("W Concentration", fontsize=fontsize)
	ax.set_ylim([0, 1e-5])
elif fom == "avg_nw":
	ax.set_ylabel("W Density (m-3)", fontsize=fontsize)

ax.legend()
fig.tight_layout()
fig.show()

# Used below
def eck_yield_ne_on_w(temp, charge):
	
	# Fitting parameters
	lamb = 0.828
	q = 2.552
	mu = 1.9534
	Eth = 38.6389
	eps = 1.19e5
	Esb = 8.68
	Esb_gam = 24.35

	# Constants and inputs
	ab = 0.0529177
	Z1 = 10
	M1 = 20.18
	Z2 = 74
	M2 = 183.84
	
	# Calculated quantities
	Eimp = 5 * temp * charge
	aL = (9 * np.pi**2 / 128)**(1/3) * ab * (Z1**(2/3) + Z2**(2/3))**(-1/2)
	El = Esb / eps
	Ered = Eimp * (M2 / (M1 + M2)) * aL / (Z1 * Z2 * 1.44)
	S_n_KrC = 0.5 * np.log(1 + 1.2288 * Ered) / (Ered + 0.1728 * np.sqrt(Ered) 
		+ 0.008 * Ered**0.1504)

	# Finally, the yield
	Y = q * S_n_KrC * (Eimp / Eth - 1)**mu / (lamb + (Eimp / Eth - 1)**mu)
	return Y

# Now we want to incorporate the linear independent variables for the figure
# of merit. Some notes:
#  - The yield in the sims was 0.001, this corresponds to a is a T=12.85 eV Ne+
#    ion entering the sheath (Eimp = 64.25 eV) at 1% the D flux.
#  - We can calculate the yield, Y2,  for a given Ne temp/charge using the above
#    function. The DIVIMP results would then scale by Y2 / 0.001, assuming
#    the concentration of Ne is the same. 
#  - The W density will just scale linearly with Ne concentration
# So with that in mind, our dependent variables are:
Y0 = 0.001  # 12.85 eV Ne+ on W
c0 = 0.01
T_Ne = np.linspace(5, 50, 10)
Z_Ne = np.arange(1, 14)
c_Ne = np.linspace(0.001, 0.1, 10)

# Then we scan through each, deriving what the total scaling factor is with
# each possible combination.
scalars = []
fom_vals_sc = []
fom_stds_sc = []
xvalues = []
for i in range(0, len(fom_vals)):
	for T in T_Ne:
		for Z in Z_Ne:
			
			# Calculate yield
			Y = eck_yield_ne_on_w(T, Z)

			for c in c_Ne:

				# Calculate scalar, append to list
				scalar = (Y / Y0) * (c / c0)
				scalars.append(scalar)
				fom_vals_sc.append(fom_vals[i] * scalar)

				# Build array that will be passed into curve fitting below
				xvalues.append([dperp[i], T, Z, c])

# Easier to work with numpy arrays
scalars = np.array(scalars)
fom_vals_sc = np.array(fom_vals_sc)
fom_stds_sc = np.array(fom_stds_sc)
xvalues = np.array(xvalues)

# Clean data, removing nan's.
mask = ~np.isnan(fom_vals_sc)
xvalues = xvalues[mask]
fom_vals_sc = fom_vals_sc[mask]

# Now that we have a large enough set of values, let's perform a multivariable
# curve fit to derive a scaling law
def model(x, a, b, c):
	
	# These are: D, T, Z and c
	print(x)
	x1, x2, x3, x4 = x

	# Concentration is by definition 1.0 exponent here, so effectively remove it
	# as a free parameter.
	return x1**a * x2**b * x3**c * x4**1.0

popt, pcov = curve_fit(model, xvalues.T, fom_vals_sc,
	bounds=((-10, -10, -10), (10, 10, 10)))
fom_fit_sc = model(xvalues.T, *popt)

# Make a ylabel from the fit
xlabel = r"$D^{:.2f}$".format(popt[0])

# Still growing this plot, but it's the figure of merit plot
fig, ax = plt.subplots()
#ax.errorbar(dperps, fom_vals, fom_stds, lw=0, marker="o", markersize=5,
#	markerfacecolor="tab:red", markeredgecolor="k", elinewidth=3, ecolor="k")
ax.scatter(fom_vals_sc, fom_fit_sc)
ax.set_xlabel(xlabel)
ax.grid(alpha=0.3)
ax.set_ylabel(fom)
fig.tight_layout()
fig.show()

