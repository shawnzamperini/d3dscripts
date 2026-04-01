# Point of this script is to compare a Flan simulation in simple helical to
# a CP inserted during those conditions (A2).
import flan_plots
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# Load A2 data
a2_path = "/home/zamp/data/A2.xlsx"
a2 = pd.read_excel(a2_path)
a2u_r = a2["R-Rsep U (cm)"].values / 100
a2d_r = a2["R-Rsep D (cm)"].values / 100
a2u_w = a2["W Areal Density U (1e15 W/cm2)"].values
a2d_w = a2["W Areal Density D (1e15 W/cm2)"].values

# Load Flan data
flan_path = "/home/zamp/flandir/sh_v7/sh_v7_w_outward_full_z.nc"
with flan_plots.FlanPlots(flan_path) as fp:
	fp_x = fp.nc["geometry"]["x"][:]
	fp_z = fp.nc["geometry"]["z"][:]

	# Restrict to t range, average along y, normalize. fp_w at this point needs
	# to be index at a specific z location to plot it against x.
	tmin = 349; tmax = 399
	fp_w = fp.nc["output"]["nz"][tmin:tmax+1, :, :, :].mean(axis=0).mean(axis=1)
	fp_w /= fp_w.max()

	# Procedure to calculate vr from 2025/10/calc_flan_dperp_v2.py
	X = fp.nc["geometry"]["X"][:, :, :].data
	Y = fp.nc["geometry"]["Y"][:, :, :].data
	vX = fp.nc["output"]["v_X"][tmin:tmax+1, :, :, :].data
	vY = fp.nc["output"]["v_Y"][tmin:tmax+1, :, :, :].data

	# Scalar projection
	vr = np.zeros(vX.shape)
	er = np.zeros(vX.shape)
	for i in range(vX.shape[0]):
		scalar_proj = X * vX[i] + Y * vY[i]
		dir_scalar = np.ones(X.shape)
		dir_scalar[scalar_proj < 0] = -1
		vr[i] = dir_scalar * np.sqrt(np.square(vX[i]) + np.square(vY[i]))

	counts = fp.nc["output"]["Nz"][tmin:tmax+1, :, :, :].data
	mask = counts > 0
	vr[~mask] = np.nan

	# Average over time, ignoring zero count regions, and then over y
	vr_avg = np.nanmean(vr, axis=0)
	vr_avg_x = np.nanmean(vr_avg, axis=1)

# Exponential fit function
def exp_fit(x, y):
	"""
	Fit y = A * exp(B*x) using linear regression on log(y).
	Returns A, B.
	"""

	# Ensure arrays
	x = np.asarray(x)
	y = np.asarray(y)

	# Filter out non-positive y values (log undefined)
	mask = y > 0
	x_fit = x[mask]
	y_fit = y[mask]

	# Linear least squares on log(y)
	logy = np.log(y_fit)
	B, logA = np.polyfit(x_fit, logy, 1)

	A = np.exp(logA)
	return A, B

# x range to do a fit
xmin = 2.31; xmax = 2.35
xmask = np.logical_and(fp_x >= xmin, fp_x <= xmax)
fp_x_fitvals = fp_x[xmask]

# Find a fit value for each z location
fp_lamb_w = []
fit_a = []; fit_b = []
nz_fit = []
dndr = []
xshift = fp_x_fitvals.min()
for zidx in range(0, fp_w.shape[1]):
	fp_w_fitvals = fp_w[xmask, zidx]
	a, b = exp_fit(fp_x_fitvals-xshift, fp_w_fitvals)
	#fp_w_fit = a * np.exp(b * fp_x_fitvals)
	lamb_w = -1 / b
	print("{} Lambda W = {:.2f} cm".format(zidx, lamb_w * 100))
	fp_lamb_w.append(lamb_w * 100)  # m to cm
	fit_a.append(a)
	fit_b.append(b)

	nz_fit.append(a * np.exp(b * (fp_x_fitvals-xshift)))

	# The fitted density gradient, exponential so just multiply by b yay
	dndr.append(b * a * np.exp(b * (fp_x_fitvals-xshift)))

# Exponential fits to the CP data in the applicable regions
a2u_mask = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
a2u_a, a2u_b = exp_fit(a2u_r[a2u_mask], a2u_w[a2u_mask])
a2d_mask = [15, 16, 17, 18]
a2d_a, a2d_b = exp_fit(a2d_r[a2d_mask], a2d_w[a2d_mask])

# -----------------
# Plotting functions
# -----------------
rsep = 2.259
fontsize = 16
fig, (ax3, ax1, ax2) = plt.subplots(1, 3, figsize=(15, 5.2))

# Plot CP data along with exponential fits
ax3.scatter(a2d_r[5:-1]*100, a2d_w[5:-1], label="ITF", s=75, color="tab:purple", 
	edgecolors="k", zorder=50)
ax3.scatter(a2u_r[5:-3]*100, a2u_w[5:-3], label="OTF", s=75, color="tab:red", 
	edgecolors="k", zorder=50)
ax3.plot(a2d_r[a2d_mask]*100, a2d_a * np.exp(a2d_b * a2d_r[a2d_mask]), 
	linestyle="--", color="tab:purple", lw=3)
ax3.plot(a2u_r[a2u_mask]*100, a2u_a * np.exp(a2u_b * a2u_r[a2u_mask]), 
	linestyle="--", color="tab:red", lw=3)
ax3.legend(fontsize=fontsize, loc="upper right")
#ax3.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
ax3.set_xlabel(r"$\mathdefault{R-R_{sep}\ (cm)}$", fontsize=fontsize)
ax3.set_ylabel(r"W Areal Density ($\mathdefault{10^{15}}$ W/$\mathdefault{cm^2}$)", fontsize=fontsize)
ax3.tick_params(axis='both', which='major', labelsize=fontsize-2)
ax3.set_ylim(0, 0.6)
ax3.set_title("Fits to Radial Profiles from CP", fontsize=fontsize)

# Show the equation for the exponential fit
ax3.text(0.68, 0.68, r"$\mathdefault{y \sim e^{-x/\lambda_{W}}}$", 
	transform=ax3.transAxes, fontsize=fontsize, color="k", 
	bbox=dict(edgecolor="k", facecolor="w"))

# Show the values of lambda for each line
ax3.text(0.2, 0.7, fr"$\mathdefault{{\lambda_{{W}}}} = {-1/a2d_b*100:.2f}$ cm", 
	transform=ax3.transAxes, fontsize=fontsize, color="tab:purple")
ax3.text(0.5, 0.15, fr"$\mathdefault{{\lambda_{{W}}}} = {-1/a2u_b*100:.2f}$ cm", 
	transform=ax3.transAxes, fontsize=fontsize, color="tab:red")

# Arrow to call out shadowed data
ax3.annotate("", xy=(8.9, 0.3), xytext=(13, 0.3),
    arrowprops=dict(arrowstyle="<-", color="tab:purple", linewidth=3, 
	mutation_scale=20))
ax3.text(0.4, 0.4, "Shadowed by wall", 
	transform=ax3.transAxes, fontsize=fontsize, color="tab:purple")

colors = ["tab:pink", "tab:cyan"]
labels = ["Near target", "Between targets"]
i = 0
for zidx in [0, 9]:
	ax1.plot((fp_x-rsep)*100, fp_w[:, zidx], color="k", lw=4)
	ax1.plot((fp_x-rsep)*100, fp_w[:, zidx], label=labels[i], 
		color=colors[i], lw=3)
	ax1.plot((fp_x_fitvals-rsep)*100, fit_a[zidx] * np.exp(fit_b[zidx] * (fp_x_fitvals-xshift)),
		color=colors[i], linestyle="--", lw=3)
	i += 1
ax1.legend(fontsize=fontsize, loc="upper right")
ax1.set_xlabel(r"$\mathdefault{R-R_{sep}\ (cm)}$", fontsize=fontsize)
ax1.set_ylabel("W Density (arb.)", fontsize=fontsize)
ax1.tick_params(axis='both', which='major', labelsize=fontsize-2)
ax1.set_title("Fits to Radial Profiles from Flan", fontsize=fontsize)

# Show the values of lambda for each line
ax1.text(0.5, 0.6, fr"$\mathdefault{{\lambda_{{W}}}} = {-1/fit_b[0]*100:.2f}$ cm", 
	transform=ax1.transAxes, fontsize=fontsize, color=colors[0])
ax1.text(0.2, 0.2, fr"$\mathdefault{{\lambda_{{W}}}} = {-1/fit_b[9]*100:.2f}$ cm", 
	transform=ax1.transAxes, fontsize=fontsize, color=colors[1])

# Show the equation for the exponential fit
#ax1.text(0.3, 0.05, r"$\mathdefault{y \sim e^{-x/\lambda_{W}}}$", 
#	transform=ax1.transAxes, fontsize=fontsize, color="k", 
#	bbox=dict(edgecolor="k", facecolor="w"))

# Range we'd roughly expect lambda_W to fall within based off the CPs
A2U_lamb_w = 1.25
A2D_lamb_w = 4.30
ax2.axhspan(-1/a2u_b*100, -1/a2d_b*100, color="tab:green", alpha=0.2)
ax2.axhline(-1/a2d_b*100, color="tab:purple", lw=3, linestyle="--")
ax2.axhline(-1/a2u_b*100, color="tab:red", lw=3, linestyle="--")
ax2.plot(fp_z+5, fp_lamb_w, lw=3, color="k")
#ax2.plot(fp_z+5, fp_lamb_w, lw=2, color="tab:green")
ax2.set_ylim(0, 5.5)
ax2.set_xlabel("Distance along field line (m)", fontsize=fontsize)
ax2.set_ylabel(r"$\mathdefault{\lambda_W\ (cm)}$", fontsize=fontsize)
ax2.tick_params(axis='both', which='major', labelsize=fontsize-2)
ax2.set_title(r"Variation of $\mathdefault{\lambda_W}$ in Flan", fontsize=fontsize)

# Plot stars corresponding to the two lines in the previous plot
ax2.scatter(fp_z[0]+5, fp_lamb_w[0], s=400, marker="*", facecolors="tab:pink", 
	edgecolors="k", zorder=50)
ax2.scatter(fp_z[9]+5, fp_lamb_w[9], s=400, marker="*", facecolors="tab:cyan", 
	edgecolors="k", zorder=50)

# Arrow to call out the range
ax2.annotate("", xy=(3, -1/a2d_b*100), xytext=(3, -1/a2u_b*100),
    arrowprops=dict(arrowstyle="<->", color="black", linewidth=3, 
	mutation_scale=20))
ax2.text(0.2, 0.85, "Range from collector\nprobe measurements", 
	transform=ax2.transAxes, fontsize=fontsize, color="k")

# a) and b) labels
ax3.text(0.014, 0.941, "a)", 
	transform=ax3.transAxes, fontsize=fontsize, color="k") 
	#bbox=dict(edgecolor="k", facecolor="w"))
ax1.text(0.014, 0.941, "b)", 
	transform=ax1.transAxes, fontsize=fontsize, color="k") 
	#bbox=dict(edgecolor="k", facecolor="w"))
ax2.text(0.014, 0.941, "c)", 
	transform=ax2.transAxes, fontsize=fontsize, color="k")
	#bbox=dict(edgecolor="k", facecolor="w"))

fig.tight_layout()
fig.show()

# Plot of radial velocities. First get zidxs of the cells around the probe.
zidxs = [6, 7, 8, 9]

# Get ranges for an error line
vr_low = [vr_avg_x[i, zidxs].min() for i in range(0, vr_avg_x.shape[0])]
vr_high = [vr_avg_x[i, zidxs].max() for i in range(0, vr_avg_x.shape[0])]

#fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 4))
fig, ax1 = plt.subplots(1, 1, figsize=(5, 4))

ax1.axhline(125, color="k", linestyle="--")
ax1.plot(fp_x-rsep, vr_avg_x[:, zidxs])
ax1.fill_between(fp_x-rsep, vr_low, vr_high)
#ax1.set_yscale("log")
ax1.set_ylim([-200, 200])

"""
for i in range(len(dndr)):
	ax2.plot(fp_x_fitvals-rsep, dndr[i])

	# xmask was defined above and is the fit region mask. vr was calculated
	# separately so it wasn't applied.
	#gamma = nz_fit[i] * vr_avg_x[:, i][xmask]
	#gamma = nz_fit[i] * vr_avg_x[:, i][xmask].mean()
	#dr = dndr[i] / gamma

	vr = vr_avg_x[:, i][xmask]
	m, c = np.polyfit(fp_x_fitvals-rsep, vr, 1)


	# With an exp fit to nz, then Dr = dndr / gamma = b / vr
	# And with a linear fit to vr, Dr = b / (m*x + c)
	#dr = fit_b[i] / vr_avg_x[:, i][xmask]
	dr = fit_b[i] / (m * (fp_x_fitvals-rsep) + c)
	ax3.plot(fp_x_fitvals-rsep, dr)
"""
fig.tight_layout()
fig.show()

