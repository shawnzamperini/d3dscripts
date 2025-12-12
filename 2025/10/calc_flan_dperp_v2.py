import numpy as np
import matplotlib.pyplot as plt
import flan_plots
from tqdm import tqdm
from scipy.interpolate import UnivariateSpline


# Load FlanPlots object

def load(path, xmin=None, xmax=None, spline_fit=False, vr_s=0, nz_s=0):

	fp = flan_plots.FlanPlots(path)

	# Load data for a range of frames where we can say it's statistically
	# settled down. Can use plot_frame_xy for this.
	fstart = 300
	#fstart = 350
	fend = 399
	#fstart = 250
	#fend = 299
	zidx = 8

	# Load first frame so we can get dimensions to set an array up.
	#X, Y, vx_avg = fp.plot_frame_xy("imp_vx", fstart, 0.0, showplot=False)
	#X, Y, nz_avg = fp.plot_frame_xy("imp_density", fstart, 0.0, showplot=False)
	#X, Y, gy_avg = fp.plot_frame_xy("imp_gyrorad", fstart, 0.0, showplot=False)

	# Coordinates
	x, y, z = fp.load_cell_centers()
	#X, Y = np.meshgrid(x, y)

	# Data arrays
	X = fp.nc["geometry"]["X"][:, :, zidx].data
	Y = fp.nc["geometry"]["Y"][:, :, zidx].data
	vX = fp.nc["output"]["v_X"][fstart:fend+1, :, :, zidx].data
	vY = fp.nc["output"]["v_Y"][fstart:fend+1, :, :, zidx].data
	nz = fp.nc["output"]["nz"][fstart:fend+1, :, :, zidx].data
	#gy = fp.nc["imp_gyrorad"][fstart:fend+1, :, :, zidx].data
	gy = np.zeros(nz.shape)
	eX = fp.nc["background"]["E_X"][fstart:fend+1, :, :, zidx].data
	eY = fp.nc["background"]["E_Y"][fstart:fend+1, :, :, zidx].data
	eZ = fp.nc["background"]["E_Z"][fstart:fend+1, :, :, zidx].data
	bX = fp.nc["background"]["B_X"][fstart:fend+1, :, :, zidx].data
	bY = fp.nc["background"]["B_Y"][fstart:fend+1, :, :, zidx].data
	bZ = fp.nc["background"]["B_Z"][fstart:fend+1, :, :, zidx].data
	ne = fp.nc["background"]["ne"][fstart:fend+1, :, :, zidx].data

	# Calculate radial values. Need the polar angle first.
	#theta = np.arctan2(vY, vX)
	#vr = vX * np.cos(theta) + vY * np.sin(theta)
	#er = eX * np.cos(theta) + eY * np.sin(theta)
	#vr = np.sqrt(vX*vX + vY*vY)
	#er = np.sqrt(eX*eX + eY*eY)

	# Scalar projection
	vr = np.zeros(vX.shape)
	er = np.zeros(vX.shape)
	for i in range(vX.shape[0]):
		scalar_proj = X * vX[i] + Y * vY[i]
		dir_scalar = np.ones(X.shape)
		dir_scalar[scalar_proj < 0] = -1
		vr[i] = dir_scalar * np.sqrt(np.square(vX[i]) + np.square(vY[i]))

	# Magnetic field magnitude
	b = np.sqrt(bX*bX + bY*bY + bZ*bZ)

	# Counts at each spot so we can ignore zero count areas (set them to nan)
	counts = fp.nc["output"]["Nz"][fstart:fend+1, :, :, zidx].data
	mask = counts > 0
	vr[~mask] = np.nan
	nz[~mask] = np.nan
	gy[~mask] = np.nan

	# Average over time, ignoring zero count regions
	vr_avg = np.nanmean(vr, axis=0)
	nz_avg = np.nanmean(nz, axis=0)
	gy_avg = np.nanmean(gy, axis=0)
	er_avg = np.nanmean(er, axis=0)
	eZ_avg = np.nanmean(eZ, axis=0)
	b_avg = np.nanmean(b, axis=0)
	ne_avg = np.nanmean(ne, axis=0)

	# Then add up the rest of the frames
	#for f in tqdm(range(fstart+1, fend+1)):
	#	X, Y, vx = fp.plot_frame_xy("imp_vx", f, 0.0, showplot=False)
	#	vx_avg += vx
	#	X, Y, nz = fp.plot_frame_xy("imp_density", f, 0.0, showplot=False)
	#	nz_avg += nz
	#	X, Y, gy = fp.plot_frame_xy("imp_gyrorad", f, 0.0, showplot=False)
	#	gy_avg += gy

	# Convert to average
	#vx_avg /= (fend - fstart)
	#nz_avg /= (fend - fstart)
	#gy_avg /= (fend - fstart)
	
	# Values for cell centers
	#x, y, z = fp.load_cell_centers()

	# Now average across the y dimension for a radial average
	iy = int(len(y) / 2)
	#vx_avg_x = vx_avg.mean(axis=1)
	#vx_avg_x = vx_avg[:, iy]
	#nz_avg_x = nz_avg.mean(axis=1)
	#nz_avg_x = nz_avg[:, iy]
	#gy_avg_x = gy_avg.mean(axis=1)
	vr_avg_x = np.nanmean(vr_avg, axis=1)
	nz_avg_x = np.nanmean(nz_avg, axis=1)
	gy_avg_x = np.nanmean(gy_avg, axis=1)
	#ex_avg_x = np.nanmean(ex_avg, axis=1)
	#ey_avg_x = np.nanmean(ey_avg, axis=1)
	er_avg_x = er_avg[:, iy]
	eZ_avg_x = eZ_avg[:, iy]
	b_avg_x = b_avg[:, iy]
	ne_avg_x = ne_avg[:, iy]

	# Overwrite with spline fit data
	if spline_fit:
		vr_spline = UnivariateSpline(x, vr_avg_x, s=vr_s, k=2)
		nz_spline = UnivariateSpline(x, nz_avg_x, s=nz_s, k=2)
		vr_avg_x = vr_spline(x)
		nz_avg_x = nz_spline(x)

	# The radial flux
	gz_avg_x = vr_avg_x * nz_avg_x

	# The radial density gradient
	dnz_dx = np.gradient(nz_avg_x, x)

	# The radial diffusion coefficient
	dr = -gz_avg_x / dnz_dx

	# Similar analysis, but for the background ions and using the ExB velocity 
	# as the radial velocity
	g_avg_x = (eZ_avg_x / b_avg_x) * ne_avg_x
	dne_dx = np.gradient(ne_avg_x, x)
	d = -g_avg_x / dne_dx

	# At each radial location, calculate the average ne across the poloidal
	# dimension, then calcuate the magnitude of the fluctuations
	eZ_avg_y = np.nanmean(eZ, axis=2)
	eZ_fluc = np.zeros(eZ.shape)
	for f in range(eZ_fluc.shape[0]):
		for ix in range(eZ_fluc.shape[1]):
			for iy in range(eZ_fluc.shape[2]):
				eZ_fluc[f, ix, iy] = (eZ[f, ix, iy] - eZ_avg_y[f, ix]) / eZ_avg_y[f, ix]
				#eZ_fluc[f, ix, iy] = eZ[f, ix, iy]

	# Mask to desired range
	if (xmin is not None and xmax is not None):
		mask = np.logical_and(x > xmin, x < xmax)
		return x[mask], vr_avg_x[mask], nz_avg_x[mask], dr[mask], \
			gy_avg_x[mask], er_avg_x[mask], eZ_avg_x[mask], d[mask], \
			eZ_fluc[:, mask, :], vr[:, mask, :]
	else:
		return x, vr_avg_x, nz_avg_x, dr, gy_avg_x, er_avg_x, eZ_avg_y, d, \
			eZ_fluc, vr
"""
# Load the radial density from 3DLIM using a case that's been in use for a bit.
import sys
sys.path.append("/home/zamp/github/utk-fusion/lim/")
import LimPlots
lim_path = "/mnt/c/Users/Shawn Zamperini/Documents/oedge_lim_files/167196-a2-tor240_44-noprobe.nc"
lp = LimPlots.LimPlots(lim_path)
lpdata = lp.plot_par_rad("nz", 21, charge="all")
mididx = np.argmin(np.abs(lpdata["X"][:, 0])) - 15
rorigin = 2.282
rad_locs = rorigin - lpdata["Y"][0]
absfac = 1e15
lim_nz = lpdata["Z"][mididx].data * absfac
lim_rzs = zip(rad_locs, np.full(len(rad_locs), -0.188))
"""
# Function to normalize data using the value at the given number
def normalize_at_r(x, y, norm_x):
	
	# Create interpolation since norm_x probably not exactly in x.
	from scipy.interpolate import interp1d
	f = interp1d(x, y)

	# Get value at x to normalize to
	norm_y = f(norm_x)

	# Normalize and return
	return x, y / norm_y

#xmin = 2.31
#xmax = 2.349
xmin = 2.305
xmax = 2.345
path = "/home/zamp/flandir/sh_v6/sh_v6_w.nc"
coll_x, coll_vx, coll_nz, coll_dr, coll_gy, ex, ey, d, eZ_fluc, coll_vx_2d \
	= load(path, xmin, xmax)
s_coll_x, s_coll_vx, s_coll_nz, s_coll_dr, coll_gy, ex, ey, d, eZ_fluc, \
	coll_vx_2d = load(path, xmin, xmax, True, 1e4, 1e-8)
path = "/home/zamp/flandir/sh_v6/sh_v6_w_nocoll.nc"
no_x, no_vx, no_nz, no_dr, no_gy, ex, ey, d, eZ_fluc, no_vx_2d \
	= load(path, xmin, xmax)
s_no_x, s_no_vx, s_no_nz, s_no_dr, no_gy, ex, ey, d, eZ_fluc, no_vx_2d \
	= load(path, xmin, xmax, True, 1e4, 1e-8)

"""
# Normalize all density data since it's largely arbitrary anyways
norm_r = 2.32
coll_x, coll_nz = normalize_at_r(coll_x, coll_nz, norm_r)
no_x, no_nz = normalize_at_r(no_x, no_nz, norm_r)
rad_locs, lim_nz = normalize_at_r(rad_locs, lim_nz, norm_r) 
"""

print("OFF: Avg. Dr = {:.3f}".format(no_dr.mean()))
print("ON:  Avg. Dr = {:.3f}".format(coll_dr.mean()))

# Averages of the fluctuation data to get a line to put on the scatter plot
#bins = np.linspace(ne_fluc.min(), ne_fluc.max(), 25)
bin_min = -4500
bin_max = 4500
bins = np.linspace(bin_min, bin_max, 25)
bin_indices = np.digitize(eZ_fluc.flatten(), bins)
eZ_fluc_avg = [eZ_fluc.flatten()[bin_indices == i].mean() for i in range(1, len(bins))]
coll_vx_2d_avg = [np.nanmean(coll_vx_2d.flatten()[bin_indices == i]) for i in range(1, len(bins))]
no_vx_2d_avg = [np.nanmean(no_vx_2d.flatten()[bin_indices == i]) for i in range(1, len(bins))]


rsep = 2.259
show_raw = True
show_fit = False
fontsize = 14
color1 = "tab:red"
color2 = "tab:purple"
lw = 4
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 6))
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4))

for ax in (ax1, ax2, ax3):
	ax.axvline(2.305-rsep, color="k", linestyle="--")
	ax.tick_params(axis='both', which='major', labelsize=fontsize-2)

ax1.axhline(0.0, color="k")
if show_raw:
	ax1.plot(no_x-rsep, no_vx, color="k", lw=lw+1, linestyle="-")
	ax1.plot(no_x-rsep, no_vx, color=color1, lw=lw, linestyle="-", label="Collisions OFF")
	ax1.plot(coll_x-rsep, coll_vx, color="k", lw=lw+1, linestyle="-")
	ax1.plot(coll_x-rsep, coll_vx, color=color2, lw=lw, linestyle="-", label="Collisions ON")
if show_fit:
	ax1.plot(s_no_x-rsep, s_no_vx, color=color1, lw=lw, label="OFF")
	ax1.plot(s_coll_x-rsep, s_coll_vx, color=color2, lw=lw, label="ON")
ax1.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{v_r\ (m/s)}$", fontsize=fontsize)
ax1.legend(loc="upper left", framealpha=1)

if show_raw:
	ax2.plot(no_x-rsep, no_nz, color="k", lw=lw+1, linestyle="-")
	ax2.plot(no_x-rsep, no_nz, color=color1, lw=lw, linestyle="-")
	ax2.plot(coll_x-rsep, coll_nz, color="k", lw=lw+1, linestyle="-")
	ax2.plot(coll_x-rsep, coll_nz, color=color2, lw=lw, linestyle="-")
if show_fit:
	ax2.plot(s_no_x-rsep, s_no_nz, color=color1, lw=lw)
	ax2.plot(s_coll_x-rsep, s_coll_nz, color=color2, lw=lw)
#if True:
#	ax2.plot(rad_locs - rsep, lim_nz, label="Experimental")
ax2.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax2.set_ylabel(r"$\mathdefault{n_W\ (arb.)}$", fontsize=fontsize)
ax2.set_yscale("log")

ax3.axhline(0.0, color="k")
if show_raw:
	ax3.plot(no_x-rsep, no_dr, color="k", lw=lw+1, linestyle="-")
	ax3.plot(no_x-rsep, no_dr, color=color1, lw=lw, label="OFF", linestyle="-")
	ax3.plot(coll_x-rsep, coll_dr, color="k", lw=lw+1, linestyle="-")
	ax3.plot(coll_x-rsep, coll_dr, color=color2, lw=lw, label="ON", linestyle="-")
if show_fit:
	ax3.plot(s_no_x-rsep, s_no_dr, color=color1, lw=lw, label="OFF")
	ax3.plot(s_coll_x-rsep, s_coll_dr, color=color2, lw=lw, label="ON")
ax3.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax3.set_ylabel(r"$\mathdefault{D_r\ (m^2/s)}$", fontsize=fontsize)
ax3.set_ylim(0.1, 300)
ax3.set_yscale("log")

"""
ax4.scatter(eZ_fluc.flatten(), no_vx_2d.flatten(), s=10, alpha=0.2, color="tab:red")
ax4.scatter(eZ_fluc.flatten(), coll_vx_2d.flatten(), s=10, alpha=0.2, color="tab:purple")
ax4.axhline(0, color="k")
ax4.plot(eZ_fluc_avg, no_vx_2d_avg, color="k", lw=3)
ax4.plot(eZ_fluc_avg, no_vx_2d_avg, color="tab:red", lw=2)
ax4.plot(eZ_fluc_avg, coll_vx_2d_avg, color="k", lw=3)
ax4.plot(eZ_fluc_avg, coll_vx_2d_avg, color="tab:purple", lw=2)
ax4.set_xlim(bin_min, bin_max)
ax4.set_xlabel("EZ Fluctuation", fontsize=fontsize)
ax4.set_ylabel("vr (m/s)", fontsize=fontsize)
"""

fig.tight_layout()
fig.show()

# Separate plot comparing gyroradius to stuff
#fig, ax1 = plt.subplots()
#ax1.scatter(coll_gy, coll_dr, color=color2, label="ON")
#ax1.scatter(no_gy, no_dr, color=color1, label="OFF")
#fig.tight_layout()
#fig.show()

fontsize=16
fig, ax3 = plt.subplots(figsize=(6, 5))
ax3.axvline(2.31-rsep, color="k", linestyle="--")
ax3.tick_params(axis='both', which='major', labelsize=fontsize-2)
ax3.axhline(0.0, color="k")
if show_raw:
	#ax3.plot(no_x-rsep, no_dr, color="k", lw=lw+1, linestyle="-")
	ax3.plot(coll_x-rsep, coll_dr, color="k", lw=lw+1, linestyle="-")
	#ax3.plot(no_x-rsep, no_dr, color=color1, lw=lw, label="OFF", linestyle="-")
	ax3.plot(coll_x-rsep, coll_dr, color=color2, lw=lw, label="ON", linestyle="-")
if False:
	ax3.plot(s_no_x-rsep, s_no_dr, color=color1, lw=lw, label="OFF")
	ax3.plot(s_coll_x-rsep, s_coll_dr, color=color2, lw=lw, label="ON")
ax3.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax3.set_ylabel(r"$\mathdefault{D_r\ (m^2/s)}$", fontsize=fontsize)
ax3.set_ylim(-10, 120)
fig.tight_layout()
fig.show()
