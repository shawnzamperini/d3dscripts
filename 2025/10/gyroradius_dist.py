import flan_plots
import matplotlib.pyplot as plt
import numpy as np


# Constants
amu_to_kg = 1.66e-27
elec = -1.602e-19

# Load a case, doesn't matter which because we just want the background
path = "/home/zamp/flandir/sh_v6/sh_v6_w_inward.nc"
fp = flan_plots.FlanPlots(path)


def get_plasma_dist_over_gyro(rho, mz, charge, center_x, center_y, center_z=0.01, tidx=0):

	# Load EZ
	EZ = fp.nc["E_Z"][:].data
	t = fp.nc["t"][:].data
	x = fp.nc["x"][:].data
	y = fp.nc["y"][:].data
	z = fp.nc["z"][:].data

	# Particle located at x, y, z will gyro-orbit perpendicular to z. x, y are
	# fortunately already in meters and the gyro-orbit plane. So at x, y, 
	# get every data point along a circle of the indicated orbit radius.
	center_xidx = np.argmin(np.abs(center_x - x))
	center_yidx = np.argmin(np.abs(center_y - y))
	center_zidx = np.argmin(np.abs(center_z - z))

	# Magnetic field variables, constant so just index at t=0.
	BX = self.nc["B_X"][0, center_xidx, center_yidx, center_xidx].data
	BY = self.nc["B_Y"][0, center_xidx, center_yidx, center_xidx].data
	BZ = self.nc["B_Z"][0, center_xidx, center_yidx, center_xidx].data
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
		orbit_t = theta / omega_c_ang
		orbit_tidx = np.argmin(np.abs(orbit_t - t))

		# Get EZ value
		orbit_EZs.append(EZ[orbit_tidx, orbit_xidx, orbit_yidx, center_zidx])
	
	return orbit_EZs

imps_mz = {"He":4.00, "Li":6.94, "B":10.81, "C":12.01, "Ne":20.18, "Fe":55.85, 
	"Mo":95.95, "W":183.84}

rhos = np.linspace(0.001, 0.025, 100)
mean_orbit_EZs = []
for rho in rhos:
	orbit_EZs = get_plasma_dist_over_gyro(rho, 2.33, 0)
	mean_orbit_EZs.append(np.mean(orbit_EZs))

fig, ax1 = plt.subplots()
ax1.plot(rhos, mean_orbit_EZs)
fig.tight_layout()
fig.show()


