import flan_plots
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


fstart = 199
fend = 399

# Load for the background data
def calc(path):
	fp = flan_plots.FlanPlots(path)

	X, Y, gyrorad_avg = fp.plot_frame_xy("imp_gyrorad", fstart, 0.0, showplot=False)

	# Then add up the rest of the frames
	for f in tqdm(range(fstart+1, fend+1)):
		X, Y, gyrorad = fp.plot_frame_xy("imp_gyrorad", f, 0.0, showplot=False)
		gyrorad_avg += gyrorad

	# Convert to average
	gyrorad_avg /= (fend - fstart)

	x, y, z = fp.load_cell_centers()

	# Now average across the y dimension for a radial average
	iy = int(len(y) / 2)
	gyrorad_avg_x = gyrorad_avg.mean(axis=1)
	#gyrorad_avg_x = gyrorad_avg[:, iy]

	return x, gyrorad_avg, gyrorad_avg_x

nocoll_x, nocoll_gyro, nocoll_gyro_x = calc("/home/zamp/flandir/testcase01/gyrorad_nocoll.nc")
coll_x, coll_gyro, coll_gyro_x = calc("/home/zamp/flandir/testcase01/gyrorad_coll.nc")

# Can I estimate the width of turbulent structures?
fp_nocoll = flan_plots.FlanPlots("/home/zamp/flandir/testcase01/gyrorad_nocoll.nc")
fp_coll = flan_plots.FlanPlots("/home/zamp/flandir/testcase01/gyrorad_coll.nc")
x, y, z = fp_nocoll.load_cell_centers()
midy = int(len(y) / 2)
midz = int(len(z) / 2)
ex = fp_nocoll.nc["elec_x"][:, :, :, midz]
ey = fp_nocoll.nc["elec_y"][:, :, :, midz]

# For each cell in space time, calculate the average ex, ey value contained 
# within a circle of the gyro radius.
tidxs = [0, 100, 200, 300]
nocoll_avg_gyro_ex = np.zeros(ex.shape)
coll_avg_gyro_ex = np.zeros(ex.shape)
nocoll_avg_gyro_ey = np.zeros(ex.shape)
coll_avg_gyro_ey = np.zeros(ex.shape)
#for tidx in tqdm(range(ex.shape[0])):
for tidx in tqdm(tidxs):
	for xidx in range(ex.shape[1]):
		print(xidx)
		for yidx in range(ex.shape[2]):

			# Get gyro radius
			gyro_nocoll = fp_nocoll.nc["imp_gyrorad"][tidx, xidx, yidx, midz]
			gyro_coll = fp_coll.nc["imp_gyrorad"][tidx, xidx, yidx, midz]

			# 2D array of distances from cell coordinates
			dist = np.zeros((ex.shape[1], ex.shape[2]))
			for xidx2 in range(ex.shape[1]):
				for yidx2 in range(ex.shape[2]):
					dist[xidx2, yidx2] = np.sqrt(np.square(x[xidx] - x[xidx2]) 
						+ np.square(y[yidx] - y[yidx2]))

			# See which indices are within +/- a gyroradius
			nocoll_mask = dist < gyro_nocoll
			coll_mask = dist < gyro_coll

			# Average values within the gyroradius
			nocoll_avg_gyro_ex[tidx, xidx, yidx] = ex[tidx][nocoll_mask].mean()
			coll_avg_gyro_ex[tidx, xidx, yidx] = ex[tidx][coll_mask].mean()
			nocoll_avg_gyro_ey[tidx, xidx, yidx] = ey[tidx][nocoll_mask].mean()
			coll_avg_gyro_ey[tidx, xidx, yidx] = ey[tidx][coll_mask].mean()


rsep = 2.259
fontsize = 14
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(11, 4))

# The gyroradii
ax1.plot(nocoll_x-rsep, nocoll_gyro[:,0], color="tab:red", alpha=0.4, label="OFF")
ax1.plot(coll_x-rsep, coll_gyro[:,0], color="tab:purple", alpha=0.4, label="ON")
ax1.plot(nocoll_x-rsep, nocoll_gyro[:,1:], color="tab:red", alpha=0.4)
ax1.plot(coll_x-rsep, coll_gyro[:,1:], color="tab:purple", alpha=0.4)

# The average plasma parameters within a gyrocircle
ax2.axhline(0, color="k")
ax3.axhline(0, color="k")
for tidx in tidxs:
	ax2.plot(nocoll_x-rsep, nocoll_avg_gyro_ex[tidx,:,midy], color="tab:red")
	ax3.plot(nocoll_x-rsep, nocoll_avg_gyro_ey[tidx,:,midy], color="tab:red")
	ax2.plot(coll_x-rsep, coll_avg_gyro_ex[tidx,:,midy], color="tab:purple")
	ax3.plot(coll_x-rsep, coll_avg_gyro_ey[tidx,:,midy], color="tab:purple")

ax1.legend()
ax1.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax2.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{\rho_W\ (m)}$", fontsize=fontsize)
ax2.set_ylabel("avg Ex", fontsize=fontsize)
ax3.set_ylabel("avg Ey", fontsize=fontsize)
fig.tight_layout()
fig.show()
