# This script is adapted from one of Tess's jupyter notebooks. As of 5/5/26, 
# this needed to be on the add-neut-diag branch of Antoine's 
# personal_gkyl_scripts repository. 

import sys

# Use pygkyl from the actual repo location since we're actively modifying it
# and I don't feel like reinstalling it each time I do
sys.path.append("/global/homes/z/zamp/personal_gkyl_scripts/pygkyl")
import pygkyl
import numpy as np

# Which data to plot
plot_param = "flan_nz"
#plot_param = "flan_ne"
polproj_plot = True
polproj_movie = False

# Load Tess's Gkeyll simulation of this WEST SOL only case
simdir = "/global/cfs/projectdirs/m3739/gkeyll/gkyl_for_flan/west-sol-62104"
fileprefix='gk_west_lsn_sol_3x2v_p1'
simulation = pygkyl.simulation_configs.import_config('west', simdir, fileprefix)

# Time in micro-seconds
simulation.normalization.set('t','mus') 

# radial coordinate normalized by the minor radius (rho=r/a)
#simulation.normalization.set('x','minor radius') 
# binormal in term of reference sound Larmor radius
simulation.normalization.set('y','Larmor radius') 
# parallel angle divided by pi
simulation.normalization.set('z','pi') 
# fluid velocity moments are normalized by the thermal velocity
simulation.normalization.set('fluid velocities','thermal velocity') 
# temperatures in electron Volt
simulation.normalization.set('temperatures','eV') 
# pressures in Pascal
simulation.normalization.set('pressures','Pa') 
# energies in mega Joules
simulation.normalization.set('energies','MJ') 

# you can check the available frames for each data type like ion_M0, 
# ion_BiMaxwellian, etc.)
sim_frames = simulation.available_frames['ion_BiMaxwellianMoments'] 
print("%g time frames available (%g to %g)" \
	%(len(sim_frames),sim_frames[0],sim_frames[-1]))

# Load in Flan data
#path = "/pscratch/sd/z/zamp/flandir/west_farsol_v1/nearsol_w_source.nc"
#path = "/pscratch/sd/z/zamp/flandir/west_nearsol_v1/ot_w_source.nc"
path = "/pscratch/sd/z/zamp/flandir/west_nearsol_point/ot_point_source.nc"
print(path)
simulation.set_flandata(path)
timeframe = simulation.flanframes[-1]  # last frame
#timeframe = simulation.flanframes[0]  # first frame
#movie_frames = np.array(sim_frames[-10:])
movie_frames = [i for i in range(1, 21)]


# Settings for each plot options
plot_settings = {
	"flan_nz": {"clim":[1e-8, 1e-4], "colorScale":"log", "colorMap":"inferno"},
	"flan_ne": {"clim":[1e17, 1e20], "colorScale":"log", "colorMap":"inferno"},
	"flan_E_X": {"clim":[-1000, 1000], "colorScale":"linear", "colorMap":"coolwarm"},
	"flan_E_Y": {"clim":[-5000, 5000], "colorScale":"linear", "colorMap":"coolwarm"},
	"flan_B_X": {"clim":[0, 3], "colorScale":"linear", "colorMap":"inferno"},
	"flan_B_R": {"clim":[0, 3], "colorScale":"linear", "colorMap":"inferno"},
	"flan_v_rad": {"clim":[-1, 1], "colorScale":"linear", "colorMap":"coolwarm"},
	"flan_v_pol": {"clim":[-1, 1], "colorScale":"linear", "colorMap":"coolwarm"},
	"ne": {"clim":[1e17, 1e20], "colorScale":"log", "colorMap":"inferno"}
}

# The nodes file is actually from a higer resolution test. This makes for a
# smoother plot.
polproj = pygkyl.PoloidalProjection()
polproj.setup(simulation,nzInterp=128, nodefilename
	='//global/cfs/cdirs/m3739/gkeyll/gkyl_for_flan/west-sol-62104/Nz128/gk_west_lsn_sol_3x2v_p1-nodes.gkyl')
#polproj.setup(simulation,nzInterp=128,nodefilename=
#              '/global/cfs/projectdirs/m3739/gkeyll/gkyl_for_flan/west-sol-62104/gk_west_lsn_sol_3x2v_p1-nodes.gkyl')
polproj.set_inset(index=0, lowerCornerRelPos=(0.1,0.2), xlim=[2.1,2.2], ylim=[-.6,-0.5], zoom=2.0, markLoc=[3,4])
polproj.set_inset(index=1, lowerCornerRelPos=(0.4,0.2), xlim=[2.35,2.55], ylim=[-.6,-0.55], zoom=3.0, markLoc=[3,4])
polproj.set_inset(index=2, lowerCornerRelPos=(0.2,0.6), xlim=[2.2,2.5], ylim=[0.55,0.67], zoom=2.0, markLoc=[1,2])
xlim=[2.0,3.0]
xlim=[2.0,3.0]
ylim=[-0.7, 0.7]
clim=plot_settings[plot_param]["clim"]

#  ---- Begin toroidal angle selection ----

# --- Choose a toroidal angle phi0 (in radians) ---
#phi0 = 0.0          # example: toroidal angle = 0
phi0 = np.pi    # example: 90 degrees toroidally

# --- Choose the z location where particles start ---
# In your case, z_min is the lower sheath boundary
#z_grid = simulation.gkyl_grid['z']   # 1D array of z cell centers
z_grid = polproj.gridsN[2]
z0 = z_grid[0]                  # pick the first z cell (near z_min)

# --- Safety factor on this flux surface ---
q0 = 8.0   # your estimate; could also load from equilibrium

# --- Compute the corresponding y value ---
y0 = z0 - q0 * phi0

# --- Convert to nearest y-index ---
#y_grid = simulation.gkyl_grid['y']   # 1D array of y cell centers
y_grid = polproj.gridsN[1]
iy0 = np.argmin(np.abs(y_grid - y0))

print("Chosen toroidal angle φ0 =", phi0)
print("Corresponding y-index =", iy0, "out of", len(y_grid))
print("y value =", y_grid[iy0])

#flan_nz = simulation.get_flanfield('flan_nz', timeframe)

# Keep only the chosen y-index
#flan_nz_filtered = np.zeros_like(flan_nz)
#flan_nz_filtered[:, iy0:iy0+1, :] = flan_nz[:, iy0:iy0+1, :]

#simulation.flanfields['flan_nz'][timeframe] = flan_nz_filtered

#  ---- End toroidal angle selection ----

climInset=clim
if polproj_plot:
	polproj.plot('ne',timeFrame=sim_frames[-1],xlim=xlim,ylim=ylim,
		colorMap='inferno', clim=[1e17, 1e19], climInset=clim, show_limiter=False)
	#polproj.plot(plot_param, timeFrame=timeframe, xlim=xlim, ylim=ylim, 
	#	colorMap=plot_settings[plot_param]["colorMap"], 
	#	colorScale=plot_settings[plot_param]["colorScale"],
	#	clim=plot_settings[plot_param]["clim"], climInset=clim, show_limiter=False,
	#	yslice=None) # Added myself in hack

# Poloidal projection movie
if polproj_movie:
	print("Generating poloidal projection movie...")
	polproj.movie(plot_param, moviePrefix='/global/homes/z/zamp/figures/flan_ex_mov_', 
	#polproj.movie(plot_param, moviePrefix='/home/zamp/figures/flan_ex_mov_', 
		xlim=xlim, ylim=ylim,
		colorScale=plot_settings[plot_param]["colorScale"],
		clim=plot_settings[plot_param]["clim"], 
		colorMap=plot_settings[plot_param]["colorMap"], 
		timeFrames=movie_frames, climInset=clim, show_limiter=False)
