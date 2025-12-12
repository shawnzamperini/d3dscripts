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
	X = fp.nc["X"][:].data
	Y = fp.nc["Y"][:].data
	Z = fp.nc["Z"][:].data
	EX = fp.nc["E_X"][:].data
	EY = fp.nc["E_Y"][:].data
	EZ = fp.nc["E_Z"][:].data
	t = fp.nc["time"][:].data
	x = fp.nc["x"][:].data
	y = fp.nc["y"][:].data
	z = fp.nc["z"][:].data
	mz = fp.nc["mz"][:].data
	qz = fp.nc["qz"][:].data

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
	BX = fp.nc["B_X"][0, center_xidx, center_yidx, center_zidx].data
	BY = fp.nc["B_Y"][0, center_xidx, center_yidx, center_zidx].data
	BZ = fp.nc["B_Z"][0, center_xidx, center_yidx, center_zidx].data
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
	
# Loop through for each impurity
lambs = []
vrads = []
gyrorads = []
rhos = []
avg_EZs = []
skew_EZs = []
dists = {}
imps = ["he", "li", "b", "c", "fe", "mo", "w"]
#imps = ["mo"]  # For quick plot iterations
for imp in imps:

	# Load each impurity simulation
	print(imp.title())
	path = "/home/zamp/flandir/sh_v6/sh_v6_{}_inward.nc".format(imp)
	fp = flan_plots.FlanPlots(path)

	# z index where everything takes place
	zidx = 8

	# Get average radial profile during steady-state
	nz_ty_avg = fp.nc["nz"][min_frame:max_frame+1, :, :, zidx] \
		.mean(axis=0).mean(axis=1)
	
	# Fit an exponential to the same range for each. Really doing a linear
	# fit to the logarithm
	x = fp.nc["x"][:]
	ln_nz = np.log(nz_ty_avg)
	slope, intercept = np.polyfit(x, ln_nz, 1)

	# The 1/e decay length is 1/slope
	lambs.append(1 / slope)

	# Now we want the average gyroradius over this range
	gyrorad = fp.calc_gyrorad()
	gyrorads.append(gyrorad[min_frame:max_frame+1, :, :, zidx].mean())

	# Likewise for the average radial velocity over this range
	vX = fp.nc["v_X"][min_frame:max_frame+1, :, :, zidx]
	vY = fp.nc["v_Y"][min_frame:max_frame+1, :, :, zidx]
	vR = np.zeros(vX.shape)

	# Projection of velocity to see if it's positive or not
	X = fp.nc["X"][:, :, zidx]
	Y = fp.nc["Y"][:, :, zidx]
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

# Convert to cm
lambs = [l * 100 for l in lambs]
gyrorads = [g * 100 for g in gyrorads]

# Calculate PDF curves
pdfs = {}
for imp in dists:
	vr = dists[imp]["vR"]
	pdf_prob, bin_edges = np.histogram(vr, bins=100, density=True)
	pdf_vr = [bin_edges[i] + (bin_edges[i+1] - bin_edges[i]) / 2 for i in 
		range(len(bin_edges)-1)]
	pdfs[imp] = {"pdf_vr":pdf_vr, "pdf_prob":pdf_prob}

# Gyroradius vs. lambda
fontsize = 16
fig, ax1 = plt.subplots()
ax1.scatter(gyrorads, lambs, marker="*", s=350, color="tab:red", 
	edgecolors="k", linewidths=2, zorder=15)
ax1.set_xlabel("Gyroradius (cm)", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{\lambda}$ (cm)", fontsize=fontsize)
ax1.tick_params(axis="both", which="both", labelsize=fontsize-2)
ax1.set_xlim([None, 0.25])
ax1.set_ylim([0.8, 2.6])
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

# PDFs of radial velocity
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

