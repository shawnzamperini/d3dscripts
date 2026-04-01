import flan_plots
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skew


# Constants
amu_to_kg = 1.66e-27
elec = -1.602e-19

# Inputs
min_frame = 349
max_frame = 399

def get_plasma_dist_over_gyro(fp, gyrorad, center_x, center_y, 
	center_z=0.01, tidx=min_frame):

	# Load from netcdf file
	X = fp.nc["geometry"]["X"][:].data
	Y = fp.nc["geometry"]["Y"][:].data
	Z = fp.nc["geometry"]["Z"][:].data
	EX = fp.nc["background"]["E_X"][:].data
	EY = fp.nc["background"]["E_Y"][:].data
	EZ = fp.nc["background"]["E_Z"][:].data
	t = fp.nc["geometry"]["time"][:].data
	x = fp.nc["geometry"]["x"][:].data
	y = fp.nc["geometry"]["y"][:].data
	z = fp.nc["geometry"]["z"][:].data
	mz = fp.nc["input"]["imp_mass_amu"][:].data * amu_to_kg
	qz = fp.nc["output"]["qz"][:].data

	# Get ER
	ER = fp.calc_E_R()

	# Starting time
	start_t = t[tidx]

	# Particle located at x, y, z will gyro-orbit perpendicular to z. x, y are
	# fortunately already in meters and the gyro-orbit plane. So at x, y, 
	# get every data point along a circle of the indicated orbit radius.
	center_xidx = np.argmin(np.abs(center_x - x))
	center_yidx = np.argmin(np.abs(center_y - y))
	center_zidx = np.argmin(np.abs(center_z - z))

	# Load gyrorad and average charge
	rho = gyrorad[tidx, center_xidx, center_yidx, center_zidx]
	charge = qz[tidx, center_xidx, center_yidx, center_zidx]

	# Magnetic field variables, constant so just index at t=0.
	BX = fp.nc["background"]["B_X"][0, center_xidx, center_yidx, center_zidx].data
	BY = fp.nc["background"]["B_Y"][0, center_xidx, center_yidx, center_zidx].data
	BZ = fp.nc["background"]["B_Z"][0, center_xidx, center_yidx, center_zidx].data
	Bsq = np.square(BX) + np.square(BY) + np.square(BZ)

	# Cyclotron frequency using center B field
	mz_kg = mz * amu_to_kg
	omega_c = np.abs(charge * elec) * np.sqrt(Bsq) / mz_kg
	omega_c_ang = 2 * np.pi / omega_c

	thetas = np.linspace(0, 2 * np.pi, 100)
	orbit_EZs = []
	for theta in thetas:

		# Get index of orbit location
		orbit_x = x + rho * np.cos(theta)
		orbit_y = y + rho * np.sin(theta)
		orbit_xidx = np.argmin(np.abs(orbit_x - x))
		orbit_yidx = np.argmin(np.abs(orbit_y - y))

		# Time is just angle divided by angular frequency (fine since theta
		# only goes from 0-2pi).
		orbit_t = start_t + theta / omega_c_ang
		orbit_tidx = np.argmin(np.abs(orbit_t - t))

		# Get EZ or ER value
		#orbit_EZs.append(EZ[orbit_tidx, orbit_xidx, orbit_yidx, center_zidx])
		orbit_EZs.append(ER[orbit_tidx, orbit_xidx, orbit_yidx, center_zidx])
	
	return rho, orbit_EZs
	
fig, ax1 = plt.subplots()

# Atomic numbers
Z = {"he": 2, "li": 3, "b": 5, "c": 6, "si": 14, "fe": 26, "mo": 42, "xe": 54, 
	"w": 74}

# Loop through for each impurity
lambs = []
vrads = []
avg_nzs = []
gyrorads = []
rhos = []
avg_EZs = []
skew_EZs = []
Zs = []
dists = {}
imps = ["he", "li", "b", "c", "si", "fe", "mo", "xe", "w"]
#imps = ["he", "li", "b", "c", "ne", "fe", "w"]  # For quick plot iterations
for imp in imps:

	# Load each impurity simulation
	print(imp.title())
	Zs.append(Z[imp])
	path = "/home/zamp/flandir/sh_v6/sh_v6_{}_inward.nc".format(imp)
	fp = flan_plots.FlanPlots(path)

	# z index where everything takes place
	zidx = 8

	# Get average radial profile during steady-state
	nz_ty_avg = fp.nc["output"]["nz"][min_frame:max_frame+1, :, :, zidx] \
		.mean(axis=0).mean(axis=1)
	
	# Fit an exponential to the same range for each. Really doing a linear
	# fit to the logarithm
	x = fp.nc["geometry"]["x"][:]

	# Fit away from source
	fit_xmin = 2.305
	fit_xmax = 2.34
	fit_xmin_idx = np.where(x == x[x >= fit_xmin][0])[0][0]
	fit_xmax_idx = np.where(x == x[x <= fit_xmax][-1])[0][0]
	x_subset = x[fit_xmin_idx:fit_xmax_idx]
	nz_subset = nz_ty_avg[fit_xmin_idx:fit_xmax_idx]

	#ln_nz = np.log(nz_ty_avg)
	ln_nz = np.log(nz_subset)
	#slope, intercept = np.polyfit(x, ln_nz, 1)
	slope, intercept = np.polyfit(x_subset, ln_nz, 1)
	ax1.scatter(x, np.log(nz_ty_avg), label=imp)
	#ax1.plot(x, slope * x + intercept)
	ax1.plot(x_subset, slope * x_subset + intercept)

	# The 1/e decay length is 1/slope
	print("  lamb = {:}".format(1 / slope))
	lambs.append(1 / slope)

	# Average density somewhere inward
	#avg_nz_x = 2.32  # About 3 cm away from source
	#avg_nz_idx = np.argmin(np.abs(x - avg_nz_x))
	avg_nz = nz_ty_avg[26:44].mean()  # Average of 2.32 +/- 0.5 cm
	avg_nzs.append(avg_nz)

	# Now we want the average gyroradius over this range
	print("  Gyroradius...")
	#gyrorad = fp.calc_gyrorad()
	#avg_gyrorad = gyrorad[min_frame:max_frame+1, fit_xmin_idx:fit_xmax_idx, :, zidx].mean()
	avg_gyrorad = 0.0  # Bypassing for now, likely not gonna do this
	print("  gyrorad = {:} m".format(avg_gyrorad))
	gyrorads.append(avg_gyrorad)

	# Likewise for the average radial velocity over this range
	print("  Radial velocity...")
	vX = fp.nc["output"]["v_X"][min_frame:max_frame+1, fit_xmin_idx:fit_xmax_idx, :, zidx]
	vY = fp.nc["output"]["v_Y"][min_frame:max_frame+1, fit_xmin_idx:fit_xmax_idx, :, zidx]
	vR = np.zeros(vX.shape)

	# Projection of velocity to see if it's positive or not
	X = fp.nc["geometry"]["X"][fit_xmin_idx:fit_xmax_idx, :, zidx]
	Y = fp.nc["geometry"]["Y"][fit_xmin_idx:fit_xmax_idx, :, zidx]
	for i in range(vX.shape[0]):
		scalar_proj = X * vX[i] + Y * vY[i]
		dir_scalar = np.ones(X.shape)
		dir_scalar[scalar_proj < 0] = -1

		vR[i] = np.sqrt(np.square(vX[i]) + np.square(vY[i])) * dir_scalar
	
	vrads.append(vR.mean())

	# Load the corresponding log file that has the velocities at xidx=33 and
	# zidx = 8 during the simulation. This is a hardcode that I put into
	# Statistics::add_vels which was, for the record:
	#if (xidx == 33 && zidx == 8)  // R-Rsep = 6 cm
	#{
	#	#pragma omp critical
	#	{
	#		std::cout << std::fixed << std::setprecision(2);
	#		std::cout << "flag " << bkg.get_X(xidx, yidx, zidx) << " "
	#			<< bkg.get_Y(xidx, yidx, zidx) << " "
	#			<< bkg.get_Z(xidx, yidx, zidx) << " "
	#			<< vX << " " << vY << " " << vZ << '\n';
	#	}
	#}
	# It is likely this was run with fewer particles to avoid the file from
	# getting too large. The lines we pull from it has the syntax:
	# flag [X] [Y] [Z] [vX] [vY] [vZ]
	"""
	dist = {"vR":[]}
	fname = "/home/zamp/flandir/sh_v6/log_{}.txt".format(imp)
	try:
		with open(fname, "r") as f:
			for line in f:
				words = line.split()
				if words[0] == "flag":
					X = float(words[1])
					Y = float(words[2])
					#Z = float(words[3])
					vX = float(words[4])
					vY = float(words[5])
					#vZ = float(words[6])

					# If the projection of (vX, vY) onto (X, Y) is positive
					# then this is a positive radial velocity, otherwise
					# it is negative and inward facing.
					dir_scalar = 1
					proj = X * vX + Y * vY
					if (proj < 0): dir_scalar = -1	
		
					# Radial velocity
					vR = dir_scalar * np.sqrt(np.square(vX) + np.square(vY))
					dist["vR"].append(vR)

		dists[imp] = dist

	# If we only have some of the log files
	except:
		print("  {} not found".format(fname))
	
	# Load EZ values over a typical gyro-orbit and average
	rho, orbit_EZ_dist = get_plasma_dist_over_gyro(fp, gyrorad, 2.33, 0.0)
	rhos.append(rho)
	avg_EZs.append(np.mean(orbit_EZ_dist))
	skew_EZs.append(skew(orbit_EZ_dist))


# Calculate PDF curves
pdfs = {}
for imp in dists:
	vr = dists[imp]["vR"]
	pdf_prob, bin_edges = np.histogram(vr, bins=100, density=True)
	pdf_vr = [bin_edges[i] + (bin_edges[i+1] - bin_edges[i]) / 2 for i in 
		range(len(bin_edges)-1)]
	pdfs[imp] = {"pdf_vr":pdf_vr, "pdf_prob":pdf_prob}
"""

# Finish plot of lambda fits
ax1.legend()
ax1.set_xlabel("x")
ax1.set_ylabel("ln(nz)")
fig.tight_layout()
fig.show()

# Convert to cm
lambs = [l * 100 for l in lambs]
gyrorads = [g * 100 for g in gyrorads]


# Gyroradius vs. lambda
fontsize = 16
fig, ax1 = plt.subplots()
ax1.scatter(gyrorads, lambs, marker="*", s=350, color="tab:red", 
	edgecolors="k", linewidths=2, zorder=15)
ax1.set_xlabel("Gyroradius (cm)", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{\lambda}$ (cm)", fontsize=fontsize)
ax1.tick_params(axis="both", which="both", labelsize=fontsize-2)
#ax1.set_xlim([None, 0.25])
#ax1.set_ylim([0.8, 2.6])
for i, label in enumerate(imps):

	# Gets crowded in the low-Z area, need to move around
	if label == "c":
		xytext = (40, 20)
	else:
		xytext = (35, 35)

	ax1.annotate(label.title(), (gyrorads[i], lambs[i]), textcoords="offset pixels", 
		xytext=xytext, ha='center', 
		arrowprops={"arrowstyle":"-", "linewidth":2, "color":"k", "zorder":5}, 
		fontsize=fontsize-2)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()

# Gyroradius vs. radial velocity
fig, ax1 = plt.subplots()
ax1.scatter(gyrorads, vrads, marker="*", s=350, color="tab:red", 
	edgecolors="k", linewidths=2, zorder=15)
ax1.set_xlabel("Gyroradius (cm)", fontsize=fontsize)
ax1.set_ylabel("Radial Velocity (m/s)", fontsize=fontsize)
ax1.tick_params(axis="both", which="both", labelsize=fontsize-2)
ax1.set_xlim([None, 0.25])
#ax1.set_ylim([0.8, 2.4])
for i, label in enumerate(imps):
    ax1.annotate(label.title(), (gyrorads[i], vrads[i]), textcoords="offset pixels", 
		xytext=(35,35), ha='center', 
		arrowprops={"arrowstyle":"-", "linewidth":2, "color":"k", "zorder":5}, 
		fontsize=fontsize-2)
fig.tight_layout()
fig.show()

# Gyroradius vs. avg_nz about 3 cm from source
fig, ax1 = plt.subplots()
ax1.scatter(gyrorads, avg_nzs, marker="*", s=350, color="tab:red", 
	edgecolors="k", linewidths=2, zorder=15)
ax1.set_xlabel("Gyroradius (cm)", fontsize=fontsize)
ax1.set_ylabel("nz 3cm from source (m/s)", fontsize=fontsize)
ax1.tick_params(axis="both", which="both", labelsize=fontsize-2)
#ax1.set_xlim([None, 0.25])
#ax1.set_ylim([0.8, 2.4])
for i, label in enumerate(imps):
    ax1.annotate(label.title(), (gyrorads[i], avg_nzs[i]), textcoords="offset pixels", 
		xytext=(35,35), ha='center', 
		arrowprops={"arrowstyle":"-", "linewidth":2, "color":"k", "zorder":5}, 
		fontsize=fontsize-2)
fig.tight_layout()
fig.show()

# Z vs. avg_nz about 3 cm from source
fig, ax1 = plt.subplots()
ax1.scatter(Zs, avg_nzs, marker="*", s=350, color="tab:red", 
	edgecolors="k", linewidths=2, zorder=15)
ax1.set_xlabel("Z", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{n_z}$ 3 cm from wall (arb.)", fontsize=fontsize)
ax1.tick_params(axis="both", which="both", labelsize=fontsize-2)
ax1.set_xlim([-7, 80])
ax1.set_ylim([2e-5, 10e-5])
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
for i, label in enumerate(imps):
    ax1.annotate(label.title(), (Zs[i], avg_nzs[i]), textcoords="offset pixels", 
		xytext=(-35,35), ha='center', 
		arrowprops={"arrowstyle":"-", "linewidth":2, "color":"k", "zorder":5}, 
		fontsize=fontsize-2)
fig.tight_layout()
fig.show()

# PDFs of radial velocity
"""
fig, ax1 = plt.subplots()
colors = ["tab:pink", "tab:cyan", "tab:brown"]
for imp, color in zip(["he", "fe", "w"], colors):
	ax1.plot(pdfs[imp]["pdf_vr"], pdfs[imp]["pdf_prob"], 
		lw=4, color="k")
	ax1.plot(pdfs[imp]["pdf_vr"], pdfs[imp]["pdf_prob"], label=imp.title(), 
		lw=3, color=color)
ax1.tick_params(axis="both", which="both", labelsize=fontsize-2)
ax1.set_xlabel("Radial Velocity (m/s)", fontsize=fontsize)
ax1.set_ylabel("Probability", fontsize=fontsize)
ax1.set_xlim([-40000, 40000])
ax1.legend(fontsize=fontsize-2)
fig.tight_layout()
fig.show()
"""
