import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.interpolate import UnivariateSpline
from scipy.signal import medfilt

# Use pygkyl from the actual repo location
import sys
sys.path.append("/home/zamp/github/personal_gkyl_scripts/pygkyl")
import pygkyl

# Don't want Times New Roman, set in pygkyl
plt.rcdefaults()

# Setup simulation object with the DIII-D NT configuration
simdir = "/home/zamp/gkyldir/NT-aug25-60x60x16res"
#simdir = "/home/zamp/gkyldir/NT-current"
fileprefix = "gk_d3d_iwl_adapt_source_3x2v_p1"
#fileprefix = "gk_d3d_negD_iwl_3x2v"
simulation = pygkyl.simulation_configs.import_config("D3D_NT", simdir, fileprefix)
simulation.normalization.set("t", "mus")  # time in microseconds

# Load in Flan data
#path = '/home/zamp/flandir/nt_v1/nt_v1_coll.nc'
#path = '/home/zamp/flandir/nt_v2/nt_v2_nocoll.nc'
#path = '/home/zamp/flandir/nt_v2/nt_v2_coll.nc'
#path = '/home/zamp/flandir/nt_v2/nt_v2_core_source.nc'
#path = '/home/zamp/flandir/nt_v3/nt_v3_nocoll.nc'
path = '/home/zamp/flandir/nt_v3/nt_v3_coll.nc'
#path = '/home/zamp/flandir/nt_v3/nt_v3_lbo_w.nc'
#path = '/home/zamp/flandir/nt_v4/nt_v4_lbo.nc'
print(path)
simulation.set_flandata(path)
timeframe = simulation.flanframes[-1]
#timeframe = simulation.flanframes[25]

# For Gkeyll data
sim_frames = simulation.available_frames['field']
#timeframe = sim_frames[-1]

# Pick plot variable here
fieldname = 'flan_nz'
clim = [1e-9, 1e-6]
colorScale = "log"
colorMap = "inferno"
"""
fieldname = 'flan_plasma_pot'
clim = [-1e4, 1e4]
"""
"""
fieldname = 'flan_gradb_X'
clim = [-0.1, 0.1]
colorScale = "linear"
colorMap = "coolwarm"
"""
"""
fieldname = 'flan_elec_X'
clim = [-1e4, 1e4]
colorScale = "log"
colorMap = "coolwarm"
"""
"""
fieldname = 'flan_B_R'
clim = [-8, 8]
colorScale = "linear"
colorMap = "coolwarm"
"""
"""
fieldname = 'flan_electron_dens'
clim = [1e17, 5e19]
colorScale = "log"
colorMap = "inferno"
"""
"""
fieldname = 'flan_ion_temp'
clim = [0.1, 1000]
colorScale = "log"
colorMap = "inferno"
"""
"""
fieldname = "flan_gca_validity_time_E"
clim = [0.01, 1000]
colorScale = "log"
colorMap = "coolwarm"
"""

# To plot Gkeyll quantities
"""
fieldname = 'ne'
clim = [1e17, 5e19]
colorScale = "log"
colorMap = "inferno"
"""
"""
fieldname = 'Bmag'
clim = [-4, 4]
colorScale = "linear"
colorMap = "coolwarm"
"""

# Poloidal projection plot
polproj = pygkyl.PoloidalProjection()
polproj.setup(simulation,nzInterp=24)
if False:
	print("Generating poloidal projection...")
	polproj.plot(fieldname,timeFrame=timeframe,colorScale=colorScale,clim=clim, 
		colorMap=colorMap, show_inset=False)
	#polproj.plot(fieldname,timeFrame=sim_frames[-1],colorScale='log')

# Poloidal projection movie
#frames = list(np.arange(-25, -1, 1))
frames = simulation.flanframes[:]
#frames = list(np.arange(500, 531))
if False:
	print("Generating poloidal projection movie...")
	polproj.movie(fieldname, moviePrefix='/home/zamp/figures/flan_ex_mov_', 
		timeFrames=frames, colorScale='log', show_inset=False, clim=clim,
		colorMap=colorMap)

# Average along y at z=0.0 showing time evolution of x profile
"""
cut_dir = 'x'
cut_coord = ['avg', 0.0]
fieldnames = [fieldname]
frames = simulation.flanframes[:]
figout = []
pygkyl.plot_utils.plot_1D_time_evolution(simulation, cut_dir, cut_coord, 
	fieldnames, frames, space_time=True, plot_type='imshow', figout=figout,
	clim=clim)
"""

# Toroidal projection plot
torproj = pygkyl.TorusProjection()
torproj.setup(simulation, Nint_polproj=24)
if False:
	print("Generating toroidal projection...")
	torproj.plot(fieldname, timeFrame=timeframe, logScale=True, clim=clim,
		colorMap="inferno")

# Toroidal projection movie
if True:
	print("Generating toroidal projection movie...")
	torproj.movie(fieldname, filePrefix='/home/zamp/figures/flan_ex_mov_', 
		timeFrames=frames, logScale=True, clim=clim, colorMap="inferno")

# Impurity values at the separatrix
if False:

	# This procedure is taken by inspecting the code in 
	# pygkyl.PoloidalProjection.plot and carrying out the necessary step to
	# have the data that gets plotted on the 2D plot. Then nearest-neighbor
	# searches for the relevant data at the separatrix coordinates for a 
	# simpler plot.

	# Import needed modules
	from pygkyl.classes import Frame, TimeSerie

	# Mimic the input variables here
	timeFrames = simulation.flanframes[-30:]
	#timeFrames = simulation.flanframes[-2:]

	# Function to load the various Flan quantities and average for the frames
	def get_mean_field_RZ(fieldName):

		# We will loop through each frame, assembling the data for each to be 
		# averaged after
		field_RZs = []
		print("Loading {}...".format(fieldName))
		for timeFrame in tqdm(timeFrames):

			# Load frame information
			with Frame(polproj.sim, fieldname=fieldName, tf=timeFrame, load=True) as field_frame:
				time = field_frame.time
				vsymbol = field_frame.vsymbol
				vunits = field_frame.vunits
				toproject = field_frame.values
				frame_info = Frame(polproj.sim, fieldname=fieldName, tf=timeFrame, load=False)
			
			field_RZs.append(polproj.project_field(toproject, frame_info).flatten())
		
		# Average across time range of frames.
		return np.mean(field_RZs, axis=0)
	
	# Load mean Flan quantities
	nW_RZ = get_mean_field_RZ("flan_nz")
	vX_RZ = get_mean_field_RZ("flan_v_X")
	vY_RZ = get_mean_field_RZ("flan_v_Y")
	vZ_RZ = get_mean_field_RZ("flan_v_Z")
	#gnzX_RZ = get_mean_field_RZ("flan_imp_gradX_nz")
	#gnzY_RZ = get_mean_field_RZ("flan_imp_gradY_nz")
	#gnzZ_RZ = get_mean_field_RZ("flan_imp_gradZ_nz")
	gnzX_RZ = np.zeros(nW_RZ.shape)
	gnzY_RZ = np.zeros(nW_RZ.shape)
	gnzZ_RZ = np.zeros(nW_RZ.shape)
	#EX_RZ = get_mean_field_RZ("flan_elec_X")
	#EY_RZ = get_mean_field_RZ("flan_elec_Y")
	#EZ_RZ = get_mean_field_RZ("flan_elec_Z")
	EX_RZ = np.zeros(nW_RZ.shape)
	EY_RZ = np.zeros(nW_RZ.shape)
	EZ_RZ = np.zeros(nW_RZ.shape)

	# Radial velocity and density gradient
	vR_RZ = np.sqrt(np.square(vX_RZ) + np.square(vY_RZ))
	gnzR_RZ = np.sqrt(np.square(gnzX_RZ) + np.square(gnzY_RZ))
	ER_RZ = np.sqrt(np.square(EX_RZ) + np.square(EY_RZ))

	# R, Z of the data
	R = polproj.RIntN.flatten()
	Z = polproj.ZIntN.flatten()

	# LCFS coordinates
	Rlcfs = polproj.Rlcfs
	Zlcfs = polproj.Zlcfs

	# Indices of nearest data value to each LCFS coordinate
	print("Computing nearest indices...")
	data_idx = [np.argmin(np.square(Rlcfs[i] - R) + np.square(Zlcfs[i] - Z))	
		for i in range(len(Rlcfs))]
	
	# The LCFS data
	nW_lcfs = nW_RZ[data_idx]
	vR_lcfs = vR_RZ[data_idx]
	vZ_lcfs = vZ_RZ[data_idx]
	gnzR_lcfs = gnzR_RZ[data_idx]
	gnzZ_lcfs = gnzZ_RZ[data_idx]
	ER_lcfs = ER_RZ[data_idx]
	EZ_lcfs = EZ_RZ[data_idx]

	# Parameterize the LCFS into distance along LCFS. The way the simulation
	# was setup, s = 0 corresponds to the IMP.
	Slcfs = np.zeros(len(Rlcfs))
	for i in range(1, len(Slcfs)):
		d = np.sqrt(np.square(Rlcfs[i-1] - Rlcfs[i]) 
			+ np.square(Zlcfs[i-1] - Zlcfs[i]))
		Slcfs[i] = Slcfs[i-1] + d
	
	# Offset so the IMP (an important region) isn't split between the ends of
	# the plot
	offset = 1.0
	smax = Slcfs.max()
	for i in range(len(Slcfs)):
		if (Slcfs[i] + offset > smax):
			Slcfs[i] = Slcfs[i] + offset - smax
		else:
			Slcfs[i] += offset
	
	# Need to sort it
	sort_idx = np.argsort(Slcfs)
	Slcfs = Slcfs[sort_idx]
	nW_lcfs = nW_lcfs[sort_idx]
	vR_lcfs = vR_lcfs[sort_idx]
	vZ_lcfs = vZ_lcfs[sort_idx]
	gnzR_lcfs = gnzR_lcfs[sort_idx]
	gnzZ_lcfs = gnzZ_lcfs[sort_idx]
	ER_lcfs = ER_lcfs[sort_idx]
	EZ_lcfs = EZ_lcfs[sort_idx]
	Rlcfs = Rlcfs[sort_idx]
	Zlcfs = Zlcfs[sort_idx]
	
	# Calculate the velocity and density gradient component perpendicular to 
	# each LCFS segment 
	vperp_lcfs = np.zeros(len(Rlcfs))
	gnzperp_lcfs = np.zeros(len(Rlcfs))
	Epol_lcfs = np.zeros(len(Rlcfs))
	for i in range(len(Rlcfs)):
		
		# Normal vector. Take advantage of simulation that Z=0 is exactly 
		# in the middle. So to get the normal vector facing inwards, we negate
		# the Z component of the LCFS segment if Z > 0, and vice-versa.

		# Components of LCFS segment
		if (i == len(Rlcfs) - 1):
			dR = Rlcfs[0] - Rlcfs[i]
			dZ = Zlcfs[0] - Zlcfs[i]
		else:
			dR = Rlcfs[i+1] - Rlcfs[i]
			dZ = Zlcfs[i+1] - Zlcfs[i]

		# Normalized normal vector, first by making Z component negative
		Nmag = np.sqrt(np.square(dR) + np.square(dZ))	
		NR = dR / Nmag
		NZ = -dZ / Nmag

		# Vector from midpoint pointing towards OUTSIDE of plasma
		centerR = 1.8
		centerZ = 0.0
		center_dR = (NR - centerR)
		center_dZ = (NZ - centerZ)

		# Midpoint of LCFS segment
		midR = dR / 2
		midZ = dZ / 2

		# If projection of normal onto vector pointing towards the center is
		# negative, then the inward facing normal will be the one from negating
		# the R component and not the Z, so we fix that if necessary here.
		scalar = (center_dR * midR + center_dZ * midZ) / \
			np.sqrt(np.square(center_dR) + np.square(center_dZ))
		if (scalar < 0):
			NR = -NR
			NZ = -NZ

		#print("{:.2f}: {:.5f} {:.5f}  N = {:.2f} {:.2f}".format(i, midR, midZ, NR, NZ))

		# Inward component of velocity
		vmag = np.sqrt(np.square(vR_lcfs[i]) + np.square(vZ_lcfs[i]))
		vperp_lcfs[i] = vR_lcfs[i] * NR + vZ_lcfs[i] * NZ

		# Same thing for density gradient, 
		gnzmag = np.sqrt(np.square(gnzR_lcfs[i]) + np.square(gnzZ_lcfs[i]))
		gnzperp_lcfs[i] = gnzR_lcfs[i] * NR + gnzZ_lcfs[i] * NZ

		# Unit vector along the LCFS segment, which defines what is called
		# the poloidal direction (technically not really poloidal, moreso
		# just parallel to the poloidal projection of the field line).
		NR_pol = dR / Nmag
		NZ_pol = dZ / Nmag

		# Projection of the electric field on the LCFS (i.e., the poloidal
		# electric field). Positive value is along the line (clockwise?)
		Epol_lcfs[i] = ER_lcfs[i] * NR_pol + EZ_lcfs[i] * NZ_pol

	# Flux entering the core
	gWin_lcfs = nW_lcfs * vperp_lcfs

	# Radial cross-separatrix diffusion coefficient
	dperp_lcfs = - gnzperp_lcfs / gWin_lcfs

	# Multiply by the toroidal area of the surface for particle rate entering
	# core. First calculate length of each flux surface. For the last element,
	# since these values are cell center values, we just approximate the final
	# L value to be the same as the previous.
	Llcfs = [Slcfs[i+1] - Slcfs[i] for i in range(len(Slcfs)-1)]
	Llcfs.append(Llcfs[-1])
	Alcfs = 2 * np.pi * Rlcfs * np.array(Llcfs)

	# Particle in rate
	phiWin_lcfs = gWin_lcfs * Alcfs

	# Normalize since no meaningful units here yet
	gWin_lcfs_norm = gWin_lcfs / np.abs(gWin_lcfs).max()
	nW_lcfs_norm = nW_lcfs / np.abs(nW_lcfs).max()

	# Smoothing
	spline = UnivariateSpline(Slcfs, medfilt(gWin_lcfs_norm, 17), s=1)
	gWin_lcfs_norm_filt = spline(Slcfs)

	simp = offset
	somp = smax / 2 + offset

	loop_data = list(zip([gWin_lcfs, nW_lcfs, phiWin_lcfs, gnzperp_lcfs, 
		dperp_lcfs, Epol_lcfs], 
		["Core W Influx (arb./m2/s)", "W Density (arb.)", "W Rate (arb./s)",
		"dnz/dperp", "Dperp", "Epol (V/m)"]))

	# Plot
	fontsize = 16

	for norm_data, label in loop_data:

		fig, ax1 = plt.subplots(figsize=(6,5))
		ax1.axhline(0.0, color="k")
		ax1.axvline(simp, color="k", linestyle="--")
		ax1.axvline(somp, color="k", linestyle="--")
		ax1.plot(Slcfs, norm_data, lw=3, color="k")
		ax1.plot(Slcfs, norm_data, lw=2, color="tab:red")
		#ax1.set_ylim([-0.5, 0.3])
		ax1.set_xlabel("Distance along LCFS (m)", fontsize=fontsize)
		ax1.set_ylabel(label, fontsize=fontsize)
		fig.tight_layout()
		fig.show()

# LBO simulation plots(s)
if False:

	# Import needed modules
	from pygkyl.classes import Frame, TimeSerie

	# Mimic the input variables here
	timeFrames = simulation.flanframes
	#timeFrames = simulation.flanframes[-10:]

	# Function to load the various Flan quantities and return a list of each frame
	def get_field_RZs(fieldName):

		# We will loop through each frame, assembling the data for each to be 
		# averaged after
		times = []
		field_RZs = []
		print("Loading {}...".format(fieldName))
		for timeFrame in tqdm(timeFrames):

			# Load frame information
			with Frame(polproj.sim, fieldname=fieldName, tf=timeFrame, load=True) as field_frame:
				time = field_frame.time
				vsymbol = field_frame.vsymbol
				vunits = field_frame.vunits
				toproject = field_frame.values
				frame_info = Frame(polproj.sim, fieldname=fieldName, tf=timeFrame, load=False)
			
			times.append(time)
			field_RZs.append(polproj.project_field(toproject, frame_info))
		
		return times, field_RZs
	
	# Load mean Flan quantities
	times, nW_RZ = get_field_RZs("flan_imp_density")

	# R, Z of the data
	R = polproj.RIntN
	Z = polproj.ZIntN

	# LCFS coordinates
	Rlcfs = polproj.Rlcfs
	Zlcfs = polproj.Zlcfs

	# Location of LBO source: LBO (X,Y,Z) = (2.11397, -0, -0.00146453)
	lbo_r = 2.11397
	lbo_z = -0.00146453

	# Got this from Copilot
	def find_cell_indices(X, Y, x, y):

		from matplotlib.path import Path

		M, N = X.shape[0] - 1, X.shape[1] - 1

		for i in range(M):
			for j in range(N):

				# Define the quadrilateral corners (counter-clockwise)
				quad = [
					(X[i, j],     Y[i, j]),
					(X[i+1, j],   Y[i+1, j]),
					(X[i+1, j+1], Y[i+1, j+1]),
					(X[i, j+1],   Y[i, j+1])
				]

				path = Path(quad)
				if path.contains_point((x, y)):
					return i, j

		print("Point not found!")
		return None  # Point not found

	# Find nearest index in R,Z data
	lbo_i, lbo_j = find_cell_indices(R, Z, lbo_r, lbo_z)
	print("(R, Z) = ({}, {})".format(R[lbo_i, lbo_j], Z[lbo_i, lbo_j]))

	# Assemble list of density data points at this location
	lbo_nW = [nW_RZ[t][lbo_i, lbo_j] for t in range(len(nW_RZ))]

	# Plot
	fig, ax = plt.subplots()
	ax.plot(times, lbo_nW)
	ax.set_xlabel("Time ()")
	ax.set_ylabel("W Density (arb.)")
	fig.tight_layout()
	fig.show()
