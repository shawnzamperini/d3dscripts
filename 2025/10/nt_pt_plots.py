import numpy as np
import matplotlib.pyplot as plt

# Use pygkyl from the actual repo location
import sys
sys.path.append("/home/zamp/github/personal_gkyl_scripts/pygkyl")
import pygkyl

# Inputs here
which = "PT"  # NT or PT
plot_param = "ne"
plot_polproj = True

# Don't want Times New Roman, set in pygkyl
plt.rcdefaults()

# Setup simulation object with the DIII-D NT/PT configuration
fileprefix = "gk_d3d_iwl_adapt_source_3x2v_p1"
if which == "PT":
	simdir = "/home/zamp/gkyldir/PT-restart-Oct25"
	simulation = pygkyl.simulation_configs.import_config("D3D_PT", simdir, 
		fileprefix)
	path = '/home/zamp/flandir/pt_v5/pt_v5.nc'
elif which == "NT":
	simdir = "/home/zamp/gkyldir/NT-restart-Oct25"
	simulation = pygkyl.simulation_configs.import_config("D3D_NT", simdir, 
		fileprefix)
	path = '/home/zamp/flandir/nt_v5/nt_v5.nc'

# Time in microseconds
simulation.normalization.set("t", "mus")

# Load in Flan data
print(path)
simulation.set_flandata(path)
timeframe = simulation.flanframes[-1]

# For Gkeyll data
sim_frames = simulation.available_frames['field']

# Settings for each plot options
plot_settings = {
	"flan_nz": {"clim":[1e-9, 1e-6], "colorScale":"log", "colorMap":"inferno"},
	"flan_E_X": {"clim":[-1000, 1000], "colorScale":"linear", "colorMap":"coolwarm"},
	"ne": {"clim":[1e16, 1e20], "colorScale":"log", "colorMap":"inferno"}
}

# Poloidal projection plot
polproj = pygkyl.PoloidalProjection()
polproj.setup(simulation,nzInterp=24)
if plot_polproj:
	print("Generating poloidal projection...")
	polproj.plot(plot_param, timeFrame=timeframe, 
		colorScale=plot_settings[plot_param]["colorScale"],
		clim=plot_settings[plot_param]["clim"], 
		colorMap=plot_settings[plot_param]["colorMap"], 
		show_inset=False)
