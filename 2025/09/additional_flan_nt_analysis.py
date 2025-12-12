import flan_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import skewnorm
from tqdm import tqdm


# Flan simulation. Below we assume that we can get a full poloidal profile
# by picking an arbitrary y value (say, the halfway point), pick a radial (x)
# coordiante, say at the separatrix, and then plot it against the poloidal
# angle (z, 0-2pi). 
ncpath = "/home/zamp/flandir/nt_v2/nt_v2_coll.nc"
fp = flan_plots.FlanPlots(ncpath)

# Radial index just inside the LCFS (127, x=0.099609375 m)
xidx = 127

# Midpoint of the y direction
yidx = int(len(fp.nc["y"][:]) / 2)

# Radial coordinate and poloidal angle
r = fp.nc["x"][:]
theta = fp.nc["z"][:]

# Electron density (t, z) and its RMS
ne_z = fp.nc["electron_dens"][:, xidx, yidx, :] / 1e18
ne_z_avg = np.mean(ne_z, axis=0)

# Density fluctuation term
ne_z_fluc = np.array([ne_t - ne_z_avg for ne_t in ne_z])

# RMS of density fluctuation
ne_z_fluc_rms = np.sqrt(np.mean(np.square(ne_z_fluc), axis=0))

#
fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.plot(theta, ne_z_fluc_rms / ne_z_avg)
ax1.set_xlabel("theta")
ax1.set_ylabel("ne_fluc_rms/ne_avg")
fig.tight_layout()
fig.show()

# Same thing, but the radial profiles of everything at various poloidal angles
ne_xs_fluc_amp = []
skews = []
locs = []
rad_skews = []
rad_locs = []
for zidx in tqdm(range(len(theta))):
	ne_x = fp.nc["electron_dens"][:, :, yidx, zidx] / 1e18
	ne_x_avg = np.mean(ne_x, axis=0)
	ne_x_fluc = np.array([ne_t - ne_x_avg for ne_t in ne_x])
	ne_x_fluc_rms = np.sqrt(np.mean(np.square(ne_x_fluc), axis=0))
	ne_x_fluc_amp = ne_x_fluc_rms / ne_x_avg

	# Append to our list
	ne_xs_fluc_amp.append(ne_x_fluc_amp)

	# Range of x values to consider for density fluctuation PDF. Indices 
	# correspond to 0.095703 (122) - 0.104297 (133) (e.g., 1 cm around the 
	# LCFS). 
	xrange = [122, 134]
	#xrange = [127, 134]
	ne_lcfs_fluc = ne_x_fluc[:, xrange[0]:xrange[1]]

	# Discretize into bins
	counts, bin_edges = np.histogram(ne_lcfs_fluc, bins=15, density=True)

	# Fit a skewed Gaussian to the data
	skew, loc, scale = skewnorm.fit(ne_lcfs_fluc.flatten())
	skews.append(skew)
	locs.append(loc)

	# Example plot just for debugging I guess
	if zidx == 16:
		fig, ax1 = plt.subplots(figsize=(5, 4))
		bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
		ax1.plot(bin_centers, counts, drawstyle='steps-mid')
		x_fit = np.linspace(bin_edges[0], bin_edges[-1], 100)
		pdf_fit = skewnorm.pdf(x_fit, skew, loc, scale)
		ax1.plot(x_fit, pdf_fit)
		fig.tight_layout()
		fig.show()
	
	# Now, find the skewness for each location along the radial coordinate
	# (we are thinking about Filip's type of plots here)
	rad_skew = []
	rad_loc = []
	for xidx in range(len(r)):	
		skew, loc, scale = skewnorm.fit(ne_x_fluc[:, xidx])
		rad_skew.append(skew)
		rad_loc.append(loc)
	rad_skews.append(rad_skew)
	rad_locs.append(rad_skew)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 8))
for zidx in range(0, len(theta), 4):
	ax1.plot(r, ne_xs_fluc_amp[zidx], label=theta[zidx])
ax1.legend()

ax2.plot(theta, locs)
ax2.set_xlabel("theta")
ax2.set_ylabel("locs")

ax3.plot(r, rad_skews[0], label=theta[0])
ax3.plot(r, rad_skews[7], label=theta[7])
ax3.plot(r, rad_skews[15], label=theta[15])
ax3.plot(r, rad_skews[23], label=theta[23])
ax3.set_xlabel("x")
ax3.set_ylabel("skewness")
ax3.set_ylim([-10, 15])
ax3.legend()

fig.tight_layout()
fig.show()
