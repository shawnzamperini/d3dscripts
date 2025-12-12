import pygkyl
import matplotlib.pyplot as plt
import numpy as np
import netCDF4


# Some options to control what plots are made
plot_movie = False
plot_sep_rad_flux = False
plot_line_density = False

# Don't want Times New Roman, set in pygkyl
plt.rcdefaults()

# Load netCDF4 file since we aren't in the flan environment
nc_path = "/home/zamp/flandir/iwl_test/iwl_test.nc"
nc = netCDF4.Dataset(nc_path)

# Setup commands
simdir = "/home/zamp/gkyldir/NT-current/"
fileprefix = "gk_d3d_negD_iwl_3x2v"
simulation = pygkyl.simulation_configs.import_config("D3D_NT", simdir, fileprefix)
polproj = pygkyl.PoloidalProjection()
polproj.setup(simulation,nzInterp=24)
sim_frames = simulation.available_frames['field']

# The plot command has been hijacked to show Flan data
R, Z, nW = polproj.plot('pi', timeFrame=[100], colorScale='log', 
	clim=[5e-11, 1e-6], inset=False)
#plt.show()

# Movie plot
if plot_movie:
	frames = list(np.arange(0, 75))
	polproj.movie('ne', moviePrefix='/home/zamp/figures/flan_ex_mov_', 
		timeFrames=frames, colorScale='log', inset=False, clim=[5e-11, 1e-6])

# Radial flux around the separatrix. Bounding x indices are 128, 129, separatrix
# is the border between these cells.
if plot_sep_rad_flux:
	ix0 = 128
	ix1 = 129
	pass

# Density along a given line
if plot_line_density:
	
	# Calculate R, Z coordinate of midpoint of each cell
	Rmids = np.zeros(R.shape)
	Zmids = np.zeros(R.shape)
	for i in range(R.shape[0]-1):
		for j in range(R.shape[1]-1):
			Rmids[i, j] = (R[i, j] + R[i+1, j] + R[i, j+1] + R[i+1, j+1]) / 4.0	
			Zmids[i, j] = (Z[i, j] + Z[i+1, j] + Z[i, j+1] + Z[i+1, j+1]) / 4.0	
	
	# Define a line of points from the start and end of measurement coords.
	line_start = (2.1, 0.0)
	line_end = (2.205, 0.0)
	ndata = 200
	line_R = np.linspace(line_start[0], line_end[0], ndata)
	line_Z = np.linspace(line_start[1], line_end[1], ndata)

	# Load each case, with or without collisions
	nc_path_nocoll = "/home/zamp/flandir/iwl_test/iwl_test_nocoll.nc"
	nc_path_coll = "/home/zamp/flandir/iwl_test/iwl_test_coll.nc"
	R, Z, nW_nocoll = polproj.plot('pi', timeFrame=[100], colorScale='log', 
		clim=[1e-9, 1e-6], inset=False, flanpath=nc_path_nocoll)
	R, Z, nW_coll = polproj.plot('pi', timeFrame=[100], colorScale='log', 
		clim=[1e-9, 1e-6], inset=False, flanpath=nc_path_coll)
	
	# Load additional frames to average
	add_frames = list(range(95, 100))
	for f in add_frames:
		R, Z, tmp_nW_nocoll = polproj.plot('pi', timeFrame=[f], colorScale='log', 
			clim=[1e-9, 1e-6], inset=False, flanpath=nc_path_nocoll)
		R, Z, tmp_nW_coll = polproj.plot('pi', timeFrame=[f], colorScale='log', 
			clim=[1e-9, 1e-6], inset=False, flanpath=nc_path_coll)
		nW_coll += tmp_nW_coll
		nW_nocoll += tmp_nW_nocoll
	nW_coll /= (len(add_frames) + 1)
	nW_nocoll /= (len(add_frames) + 1)

	# For each line coordinate, find nearest data point save data.
	data_nocoll = np.zeros(ndata)
	data_coll = np.zeros(ndata)
	for d in range(ndata):
		tmpR = line_R[d]
		tmpZ = line_Z[d]

		# Calculate distance to each point
		dist = np.sqrt(np.square(tmpR - Rmids) + np.square(tmpZ - Zmids))

		# Get data at smallest distance
		min_index = np.unravel_index(np.argmin(dist), dist.shape)
		data_nocoll[d] = nW_nocoll[min_index]
		data_coll[d] = nW_coll[min_index]
	
	# Plot at OMP
	fig, ax1 = plt.subplots(figsize=(7, 6))
	ax1.axvline(2.17, color="k", linestyle="--", lw=2)
	ax1.axvline(2.20, color="k", linestyle="-", lw=2)
	ax1.plot(line_R, data_nocoll, lw=3, color="k")
	ax1.plot(line_R, data_nocoll, lw=2, color="tab:red", label="OFF")
	ax1.plot(line_R, data_coll, lw=3, color="k")
	ax1.plot(line_R, data_coll, lw=2, color="tab:purple", label="ON")
	ax1.set_xlabel("R (m)", fontsize=14)
	ax1.set_ylabel("W Density (arb.)", fontsize=14)
	ax1.set_ylim([1e-10, 1e-5])
	ax1.set_yscale("log")
	ax1.grid(alpha=0.3)
	ax1.legend(fontsize=16)
	fig.tight_layout()
	fig.show()
