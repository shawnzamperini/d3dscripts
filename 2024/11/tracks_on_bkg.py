import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import flan_plots
from matplotlib.colors import Normalize, LogNorm
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import PillowWriter


# Various options
data_name = "elec_y"
frame_start = 0
frame_end = 49
z0 = 0.0
norm_type = "linear"
vmin = -5000
vmax = 5000
cmap = "coolwarm"
rsep = 0.0
xlabel = "R (m)"
ylabel = "Binormal (m)"
g_fontsize = 14
cbar_label = "Ey (V/m)"
save_path = "W3coll"
#save_path = None
sheet_name = "W3+ coll"
#dt = 1e-7

# We just need this for the background data. Impurities are being plotted
# a different way here.
path = "/home/zamp/flandir/testcase01/coll_on_v3.nc"
fp = flan_plots.FlanPlots(path)

# This contains the W particle tracks for different charge states.
xl_path = "/mnt/g/My Drive/Research/Documents/2024/flan_w_tracks.xlsx"
df = pd.read_excel(xl_path, sheet_name="{}".format(sheet_name))

# Assign particle ID as index, useful for grouping later
df.set_index("id", inplace=True)

# Create another frame index, where instead of the background time resolution
# it is the time resolution of the impurities (dt). 
#df["imp_frame"] = ((df["t"] - fp.nc["time"][0]) / dt).astype(int)

# Only works for a single impurity with variable time step. I've changed this
# so the whole script don't work beyond one particle now.
df["imp_frame"] = np.arange(0, len(df))

# Determine frame range from the impurity timestep perspective.
imp_frame_start = df["imp_frame"].min()
imp_frame_end = df["imp_frame"].max()

# A dictionary used to map imp_frame --> bkg_frame
imp_to_bkg = {}
for f in df["imp_frame"].unique():
	if f not in imp_to_bkg.keys():
		imp_to_bkg[f] = df[df["imp_frame"] == f]["tidx"].iloc[0]

# Create a list of smaller DataFrames that contain each particle's average
# at each time frame it encountered.
#avg_imps = []
#for i in df.index.unique():
#	avg_imps.append(df.loc[i].groupby("tidx").mean())

# Create plot that will be an animation. This is all copied from flan_plots. 
# Setup the first frame.
X, Y, data_xy = fp.plot_frame_xy(data_name, frame_start, z0, showplot=False)

# Removed long comment, but copied from plot_frame_xy. 
fig, ax1 = plt.subplots()
div = make_axes_locatable(ax1)
cax = div.append_axes("right", "5%", "5%")
norm = fp.get_norm(data_xy, norm_type, vmin=vmin, vmax=vmax)
mesh = ax1.pcolormesh(X-rsep, Y, data_xy.T, cmap=cmap, norm=norm)
cbar = fig.colorbar(mesh, cax=cax)
ax1.set_facecolor("grey")
ax1.set_aspect("equal")
ax1.set_xlabel(xlabel, fontsize=g_fontsize)
ax1.set_ylabel(ylabel, fontsize=g_fontsize)
ax1.set_title("Frame {}".format(frame_start), fontsize=g_fontsize)
if cbar_label is None:
	cbar.set_label(data_name, fontsize=g_fontsize)
else:
	cbar.set_label(cbar_label, fontsize=g_fontsize)
fig.tight_layout()

# Define animation function
def animate(i):

	# For this impurity frame, find what background frame is being used
	imp_frame = imp_frame_start + i
	bkg_frame = imp_to_bkg[imp_frame]

	# Call for the next frame. i is being passed in zero-indexed, so to
	# get the frame we offset it from the start_frame.
	X, Y, data_xy = fp.plot_frame_xy(data_name, bkg_frame, z0,
	   showplot=False)

	# Time at frame, starting it at zero
	time = fp.nc["time"][bkg_frame - frame_start] - fp.nc["time"][0]

	norm = fp.get_norm(data_xy, norm_type, vmin=vmin, vmax=vmax)
	ax1.clear()
	mesh = ax1.pcolormesh(X-rsep, Y, data_xy.T, cmap=cmap, norm=norm)
	ax1.set_title("{:.2f} us".format(time * 1e6))
	ax1.set_xlabel(xlabel, fontsize=g_fontsize)
	ax1.set_ylabel(ylabel, fontsize=g_fontsize)
	cax.cla()
	fig.colorbar(mesh, cax=cax)
	if cbar_label is None:
		cbar.set_label(data_name, fontsize=g_fontsize)
	else:
		cbar.set_label(cbar_label, fontsize=g_fontsize)

	# Plot particle position. Loop through all particles, adding a marker 
	# whenever a particle was in the current frame. The marker is at the
	# particle's average location during this frame.
	#for imp in avg_imps:
	#	try:
	#		s = imp.loc[i]
	#		ax1.scatter(s["x"], s["y"], marker=".", s=75, color="white", 
	#			edgecolors="k")
	#	except KeyError:
	#		pass

	# Plot the track so far. When the y coordinate wraps around it'll plot
	# a vertical line, we don't want that so we break this into a number of
	# different lines, broken whenever the data switches sign by more than
	# some reasonable amount so as to not care when it's crossing zero.
	tmp_df = df[df["imp_frame"] <= imp_frame]
	tmp_x = []
	tmp_y = []
	for i in range(0, len(tmp_df)-1):
		tmp_x.append(tmp_df["x"].iloc[i])
		tmp_y.append(tmp_df["y"].iloc[i])

		# If we wrap around in the next frame, start a new line
		if (np.abs(tmp_df["y"].iloc[i]) > 0.005 
			and ((tmp_df["y"].iloc[i] > 0 and tmp_df["y"].iloc[i+1] < 0)
			or (tmp_df["y"].iloc[i] < 0 and tmp_df["y"].iloc[i+1] > 0))
			or i == len(tmp_df)-2):
			ax1.plot(tmp_x, tmp_y, lw=3, color="k", zorder=5)
			tmp_x = []
			tmp_y = []

	# Plot particle(s) position where the data exists for this imp frame.
	tmp_df = df[df["imp_frame"] == imp_frame]
	for col, s in tmp_df.iterrows():
		ax1.scatter(s["x"], s["y"], marker=".", s=150, color="white", 
			edgecolors="k", zorder=15)

	return mesh,

# Generate animation
#anim = animation.FuncAnimation(fig, animate, frames=frame_end-frame_start)
anim = animation.FuncAnimation(fig, animate, 
	frames=imp_frame_end-imp_frame_start)

# If save_path is provided, save the file there as a gif.
if (save_path is not None):
	print("Saving animation...")
	anim.save(save_path + ".gif", writer=PillowWriter(fps=60))

plt.show()
