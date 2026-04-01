import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
import flan_plots


# Make plots for only the true cases
make_curv = False
make_exb  = False
make_dBdy = False
make_dEdt = False
make_ff   = True
make_num  = False

# Load Excel sheet. Need to mount G Drive:
# sudo mount -t drvfs G: /mnt/g
path = "/mnt/g/My Drive/Research/Documents/2026/exb_test.xlsx"

# --------------
# Curvature Drift
# ---------------
if make_curv:
	print("Loading curvature drifts test...")
	df = pd.read_excel(path, sheet_name="test_gyro_curv_hires")

	vpar_init = 5000
	t = df["t"].values
	x = df["x"].values
	y = df["y"].values
	z = df["z"].values

	# Normalize or whatever
	t = (t - t.min()) * 1e6
	y = y * 1000  # m to mm

	# Estimate of the analytical curvature drift
	q = 1.602e-19
	m = 183.84 * 1.66e-27
	omega_c = q / m  # Assumes B = 1
	avg_x = (x.min() + x.max()) / 2.0
	calc_curv = -vpar_init**2 / omega_c / avg_x

	# Find the local maxima so we can track the drift motion
	peaks, properties = find_peaks(y) # Extract peak times and values 
	t_peaks = t[peaks] 
	y_peaks = y[peaks]
	vk, b = np.polyfit(t_peaks, y_peaks, 1)
	tfit = np.linspace(t.min(), t.max(), 100)
	yfit = vk * tfit + b
	print("Curvature drift")
	print("  Measured:   {:.2f}".format(vk*1e3))
	print("  Calculated: {:.2f}".format(calc_curv))

	fontsize = 14
	fig, ax1 = plt.subplots(figsize=(5, 4))
	#ax1.axhline(0.0, color="k")
	ax1.plot(t, y, zorder=5, lw=3, color="k")
	ax1.plot(t, y, zorder=5, lw=2, color="tab:red")
	ax1.plot(tfit, yfit, linestyle="--", color="k", lw=2)
	ax1.scatter(t_peaks, y_peaks, zorder=15, color="tab:red", edgecolors="k", s=75)
	ax1.set_xlabel("Time (us)", fontsize=fontsize)
	ax1.set_ylabel("Z (mm)", fontsize=fontsize)
	ax1.tick_params(axis='both', which='major', labelsize=fontsize-2)
	ax1.set_ylim([-7, 10])
	ax1.text(0.33, 0.80, fr"$\mathdefault{{v_{{\kappa}} = {vk*1e3:.2f}}}$ m/s", 
		transform=ax1.transAxes, fontsize=fontsize, color="k", rotation=7.2)
	fig.tight_layout()
	fig.show()

if make_exb:
	print("Loading ExB drifts test...")
	df = pd.read_excel(path, sheet_name="test_gyro_ExB")

	vY_init = 2500
	t = df["t"].values
	x = df["x"].values
	y = df["y"].values
	z = df["z"].values

	# Normalize or whatever
	t = (t - t.min()) * 1e6
	x = (x - x.min()) * 100

	# Calculate expected ExB drift
	EY = 500
	BZ = 1.0
	calc_exb = EY / BZ

	# Find the local maxima so we can track the drift motion
	peaks, properties = find_peaks(x) # Extract peak times and values 
	t_peaks = t[peaks] 
	x_peaks = x[peaks]
	vexb, b = np.polyfit(t_peaks, x_peaks, 1)
	tfit = np.linspace(t.min(), t.max(), 100)
	xfit = vexb * tfit + b
	print("ExB drift")
	print("  Measured:   {:.2f}".format(vexb*1e4))
	print("  Calculated: {:.2f}".format(calc_exb))

	fontsize = 14
	fig, ax1 = plt.subplots(figsize=(5, 4))
	#ax1.axhline(0.0, color="k")
	ax1.plot(t, x, zorder=5, lw=3, color="k")
	ax1.plot(t, x, zorder=5, lw=2, color="tab:red")
	ax1.plot(tfit, xfit, linestyle="--", color="k", lw=2)
	ax1.scatter(t_peaks, x_peaks, zorder=15, color="tab:red", edgecolors="k", s=75)
	ax1.set_xlabel("Time (us)", fontsize=fontsize)
	ax1.set_ylabel("x (cm)", fontsize=fontsize)
	ax1.tick_params(axis='both', which='major', labelsize=fontsize-2)
	ax1.set_ylim([-0.2, 4])
	ax1.text(0.15, 0.48, fr"$\mathdefault{{v_{{ExB}} = {vexb*1e4:.2f}}}$ m/s", 
		transform=ax1.transAxes, fontsize=fontsize, color="k", rotation=25)
	fig.tight_layout()
	fig.show()

if make_dBdy:
	print("Loading Grad-B drifts test...")
	df = pd.read_excel(path, sheet_name="test_gyro_dBdy")

	vY_init = 2500
	t = df["t"].values
	x = df["x"].values
	y = df["y"].values
	z = df["z"].values

	# Normalize or whatever
	t = (t - t.min()) * 1e6
	x = (x - x.min()) * 100

	# Calculate expected ExB drift
	q = 1.602e-19
	m = 183.84 * 1.66e-27
	omega_c = q / m  # Assumes B = 1
	dBdy = 5  # 5 T/m or 0.05 T/cm
	calc_gradB = vY_init**2 / (2 * omega_c) * dBdy

	# Find the local maxima so we can track the drift motion
	peaks, properties = find_peaks(x) # Extract peak times and values 
	t_peaks = t[peaks] 
	x_peaks = x[peaks]
	vgradB, b = np.polyfit(t_peaks, x_peaks, 1)
	tfit = np.linspace(t.min(), t.max(), 100)
	xfit = vgradB * tfit + b
	print("Grad-B drift")
	print("  Measured:   {:.2f}".format(vgradB*1e4))
	print("  Calculated: {:.2f}".format(-calc_gradB))

	fontsize = 14
	fig, ax1 = plt.subplots(figsize=(5, 4))
	#ax1.axhline(0.0, color="k")
	ax1.plot(t, x, zorder=5, lw=3, color="k")
	ax1.plot(t, x, zorder=5, lw=2, color="tab:red")
	ax1.plot(tfit, xfit, linestyle="--", color="k", lw=2)
	ax1.scatter(t_peaks, x_peaks, zorder=15, color="tab:red", edgecolors="k", s=75)
	ax1.set_xlabel("Time (us)", fontsize=fontsize)
	ax1.set_ylabel("x (cm)", fontsize=fontsize)
	ax1.tick_params(axis='both', which='major', labelsize=fontsize-2)
	ax1.set_ylim([-0.1, 1.5])
	ax1.text(0.15, 0.8, fr"$\mathdefault{{v_{{\nabla B}} = {vgradB*1e4:.2f}}}$ m/s", 
		transform=ax1.transAxes, fontsize=fontsize, color="k", rotation=-8)
	fig.tight_layout()
	fig.show()

if make_dEdt:
	print("Loading polarization drifts test...")
	df = pd.read_excel(path, sheet_name="test_gyro_dEdt")

	vpar_init = 2500
	t = df["t"].values
	x = df["x"].values
	y = df["y"].values
	z = df["z"].values

	# Normalize or whatever
	t = (t - t.min()) * 1e6
	y = y * 1000  # m to mm

	# Estimate of the analytical polarization drift
	q = 1.602e-19
	m = 183.84 * 1.66e-27
	omega_c = q / m  # Assumes B = 1
	dEdt = -1000000  # -1000000 V/s or -1 V/us
	calc_pol = dEdt / (omega_c)

	# Find the local maxima so we can track the drift motion
	peaks, properties = find_peaks(y) # Extract peak times and values 
	t_peaks = t[peaks] 
	y_peaks = y[peaks]
	vpol, b = np.polyfit(t_peaks, y_peaks, 1)
	tfit = np.linspace(t.min(), t.max(), 100)
	yfit = vpol * tfit + b
	print("Polarization drift")
	print("  Measured:   {:.2f}".format(vpol*1e3))
	print("  Calculated: {:.2f}".format(calc_pol))

	fontsize = 14
	fig, ax1 = plt.subplots(figsize=(5, 4))
	#ax1.axhline(0.0, color="k")
	ax1.plot(t, y, zorder=5, lw=3, color="k")
	ax1.plot(t, y, zorder=5, lw=2, color="tab:red")
	ax1.plot(tfit, yfit, linestyle="--", color="k", lw=2)
	ax1.scatter(t_peaks, y_peaks, zorder=15, color="tab:red", edgecolors="k", s=75)
	ax1.set_xlabel("Time (us)", fontsize=fontsize)
	ax1.set_ylabel("z (mm)", fontsize=fontsize)
	ax1.tick_params(axis='both', which='major', labelsize=fontsize-2)
	ax1.set_ylim([-7, 7])
	ax1.text(0.35, 0.90, fr"$\mathdefault{{v_{{pol}} = {vpol*1e3:.2f}}}$ m/s", 
		transform=ax1.transAxes, fontsize=fontsize, color="k", rotation=-1)
	fig.tight_layout()
	fig.show()

if make_ff:
	
	# No Excel for this one, use NC file, slab geometry
	ncpath = "/home/zamp/flandir/test_gyro/test_gyro_no_EM_ff_w.nc"
	fp_w = flan_plots.FlanPlots(ncpath)
	#ncpath = "/home/zamp/flandir/test_gyro/test_gyro_BX_ff_w.nc"
	ncpath = "/home/zamp/flandir/test_gyro/test_gyro_ff_w.nc"
	fp_w_b = flan_plots.FlanPlots(ncpath)
	#ncpath = "/home/zamp/flandir/test_gyro/test_gyro_no_EM_ff_he.nc"
	#fp_he = flan_plots.FlanPlots(ncpath)

	# This parallel flow on the ions was imposed, other values included
	vi_par = 1000  # m/s in X direction, slab
	Ti = 1.0
	mi = 2.014  # amu
	ni = 1e20
	ln_alpha = 15
	Tz = Ti
	Z = 1
	#t0 = 0.0003
	t0 = 0.0
	t = fp_w.nc["geometry"]["time"][:]
	t = t - t0

	d = {}
	#for imp, fp in zip(["w", "he"], [fp_w, fp_he]):
	for imp, fp in zip(["w", "w_b"], [fp_w, fp_w_b]):

		mz = fp.nc["input"]["imp_mass_amu"][:] # amu
		mz_kg = mz * 1.66e-27 # amu

		# Average and max vZ vs. t
		vz = fp.nc["output"]["v_X"][:]
		vz[vz == 0.0] = np.nan
		vz_t_avg = np.nanmean(np.nanmean(np.nanmean(vz, axis=1), axis=1), axis=1)
		vz_t_max = fp.nc["output"]["v_X"][:].max(axis=1).max(axis=1).max(axis=1)

		# Stopping power time from Stangeby 6.35, in s
		tau_s = 1.47e13 * mz * Ti * np.sqrt(Ti / mi) / ((1 + mi / mz) * ni * Z**2 * ln_alpha)

		# Friction force
		ff_avg = mz_kg * (vi_par - vz_t_avg) / tau_s
		ff_max = mz_kg * (vi_par - vz_t_max) / tau_s

		# Simple first-order ODE to solve for vz from F=ma
		vz_calc = vi_par * (1 - np.exp(-t / tau_s))

		d[imp] = {"vz_calc": vz_calc, "vz_t_avg": vz_t_avg}

	t_us = t * 1e6

	# Mask data at end where noise blows up
	mask = t_us < 75

	fontsize = 14
	fig, ax1 = plt.subplots(figsize=(5, 4))
	ax1.plot(t_us[mask], d["w"]["vz_calc"][mask], color="k", linestyle="--", lw=2, label="Analytic")
	ax1.plot(t_us[mask], d["w"]["vz_t_avg"][mask], color="k", lw=3)
	ax1.plot(t_us[mask], d["w"]["vz_t_avg"][mask], label="No B", color="tab:red", lw=2)
	ax1.plot(t_us[mask], d["w_b"]["vz_t_avg"][mask], color="k", lw=3)
	ax1.plot(t_us[mask], d["w_b"]["vz_t_avg"][mask], label="B=1.0 T", color="tab:purple", lw=2)
	ax1.legend(fontsize=fontsize)
	ax1.set_xlabel("Time (us)", fontsize=fontsize)
	ax1.set_ylabel(r"$\mathdefault{v_X}$ (m/s)", fontsize=fontsize)
	ax1.tick_params(axis='both', which='major', labelsize=fontsize-2)

	fig.tight_layout()
	fig.show()

# -------------------
# Numerical Diffusion
# -------------------
if make_num:
	print("Loading numerical diffusion test...")
	df = pd.read_excel(path, sheet_name="test_gyro_curv_hires")

	t = df["t"].values
	x = df["x"].values
	y = df["y"].values
	z = df["z"].values
	X = df["X"].values
	Y = df["Y"].values
	Z = df["Z"].values
	R = np.sqrt(np.square(X) + np.square(Y))

	# Normalize or whatever
	t = (t - t.min()) * 1e6

	# Fit to peaks to estimate numerical diffusion
	peaks, properties = find_peaks(x) # Extract peak times and values 
	t_peaks = t[peaks] 
	x_peaks = x[peaks]
	vnum, b = np.polyfit(t_peaks, x_peaks, 1)
	tfit = np.linspace(t.min(), t.max(), 100)
	xfit = vnum * tfit + b

	# Likewise but for Cartesian
	peaks, properties = find_peaks(R, distance=10000) # Extract peak times and values 
	tR_peaks = t[peaks] 
	R_peaks = R[peaks]
	vnumR, bR = np.polyfit(tR_peaks, R_peaks, 1)
	Rfit = vnumR * tfit + bR

	fontsize = 14
	fig, ax1 = plt.subplots(figsize=(5, 4))
	ax1.plot(t, x, lw=3, color="k")
	ax1.plot(t, x, label="Field-aligned", lw=2, color="tab:red")
	ax1.plot(t, R, lw=3, color="k")
	ax1.plot(t, R, label="Cartesian", lw=2, color="tab:purple")
	ax1.plot(tfit, xfit, linestyle="--", color="k", lw=2)
	ax1.plot(tfit, Rfit, linestyle="--", color="k", lw=2)
	ax1.scatter(t_peaks, x_peaks, zorder=15, color="tab:red", edgecolors="k", s=75)
	ax1.scatter(tR_peaks, R_peaks, zorder=15, color="tab:purple", edgecolors="k", s=75)
	ax1.set_xlabel("Time (us)", fontsize=fontsize)
	ax1.set_ylabel("R (m)", fontsize=fontsize)
	ax1.legend(fontsize=fontsize)
	ax1.set_ylim([2.21, 2.34])
	ax1.tick_params(axis='both', which='major', labelsize=fontsize-2)
	ax1.text(0.35, 0.90, fr"$\mathdefault{{v_{{num}} = {vnum*1e3:.2e}}}$ m/s", 
		transform=ax1.transAxes, fontsize=fontsize, color="k", rotation=0)
	ax1.text(0.35, 0.345, fr"$\mathdefault{{v_{{num}} = {vnumR*1e3:.2f}}}$ m/s", 
		transform=ax1.transAxes, fontsize=fontsize, color="k", rotation=-33.5)

	fig.tight_layout()
	fig.show()
