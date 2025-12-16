import numpy as np
import matplotlib.pyplot as plt

# Paths are using Perlmutter locations since the NT-high-time-res files are
# kinda too large for my laptop

# Use pygkyl from the actual repo location
import sys
sys.path.append("/global/homes/z/zamp/personal_gkyl_scripts/pygkyl")
import pygkyl

# Inputs here
which = "NT"  # NT or PT
plot_param = "flan_nz"
plot_polproj = False
polproj_movie = True
movie_frames = "all"  # "all" or [min, max]

# Don't want Times New Roman, set in pygkyl
plt.rcdefaults()

# Setup simulation object with the DIII-D NT/PT configuration
fileprefix = "gk_d3d_iwl_adapt_source_3x2v_p1"
if which == "PT":
	simdir = "/global/cfs/cdirs/m3739/gkeyll/gkyl_for_flan/PT-high-time-res"
	simulation = pygkyl.simulation_configs.import_config("D3D_PT", simdir, 
		fileprefix)
	path = '/pscratch/sd/z/zamp/flandir/pt_v5/sol_source_v1.nc'
elif which == "NT":
	simdir = "/global/cfs/cdirs/m3739/gkeyll/gkyl_for_flan/NT-high-time-res"
	simulation = pygkyl.simulation_configs.import_config("D3D_NT", simdir, 
		fileprefix)
	path = '/pscratch/sd/z/zamp/flandir/nt_v5/sol_source_v1.nc'

# Time in microseconds
simulation.normalization.set("t", "mus")

# Load in Flan data
print(path)
simulation.set_flandata(path)
timeframe = simulation.flanframes[-1]

# For Gkeyll data
sim_frames = simulation.available_frames['field']

# Select frames for movies
frames = simulation.flanframes[:]
if movie_frames == "all": print("All frames in movie")
elif type(movie_frames) is list: 
	frames = frames[movie_frames[0]:movie_frames[1]]
else: 
	print("Error: Unrecognized option. movie_frames = {}".format(movie_frames))

# Settings for each plot options
plot_settings = {
	"flan_nz": {"clim":[1e-9, 1e-6], "colorScale":"log", "colorMap":"inferno"},
	"flan_E_X": {"clim":[-1000, 1000], "colorScale":"linear", "colorMap":"coolwarm"},
	"ne": {"clim":[1e17, 1e20], "colorScale":"log", "colorMap":"inferno"}
}

# Poloidal projection plot
polproj = pygkyl.PoloidalProjection()
polproj.setup(simulation, nzInterp=24)
if plot_polproj:
	print("Generating poloidal projection...")
	polproj.plot(plot_param, timeFrame=timeframe, 
		colorScale=plot_settings[plot_param]["colorScale"],
		clim=plot_settings[plot_param]["clim"], 
		colorMap=plot_settings[plot_param]["colorMap"], 
		show_inset=False)

# Poloidal projection movie
if polproj_movie:
	print("Generating poloidal projection movie...")
	polproj.movie(plot_param, moviePrefix='/global/homes/z/zamp/figures/flan_ex_mov_', 
		colorScale=plot_settings[plot_param]["colorScale"],
		clim=plot_settings[plot_param]["clim"], 
		colorMap=plot_settings[plot_param]["colorMap"], 
		timeFrames=frames, show_inset=False)
