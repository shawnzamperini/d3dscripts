import json
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import numpy as np


# Load json with data from Luca on prompt redeposition
path = "/home/zamp/documents/DIVIMP_fnp_database.json"
with open(path, "r") as f:
	data = json.load(f)

# Note: fnp is fraction of non-promptly redeposited (fnp = 1 - fprompt)

# Let's make plots. First create normalization for colormap scaling
cmap = cm.inferno
norm = LogNorm(vmin=data["W"]["ne"][0], vmax=data["W"]["ne"][-1])

# Create figure and axes, making room for a cbar
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
fig.subplots_adjust(bottom=0.15, right=0.85)

# Add a colorbar of the colormap scaling
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.90, 0.2, 0.03, 0.65])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label("ne (m-3)")

# Plot each fprompt
for s, ax in zip(["W", "Fe"], [ax1, ax2]):
	d = data[s]
	x = d["Te"]
	for i in range(0, len(d["ne"])):
		y = d["fnp"][i]
		color = cmap(norm(d["ne"][i]))
		ax.plot(x, y, color="k", lw=3)
		ax.plot(x, y, color=color, lw=2)
	ax.set_xlabel("Te (eV)")
	ax.set_ylim([0, None])

ax1.set_ylabel("1 - fprompt")
ax1.set_title("Tungsten")
ax2.set_title("Iron")

fig.show()

def print_fortran_array(arr, name, fmt="float"):

	print("real :: {}({}) = (/".format(name, len(arr)), end="")
	for i in range(len(arr)-1):
		if fmt == "float":
			print("{:.2f}, ".format(arr[i]), end="")
		elif fmt == "sci":
			print("{:.2e}, ".format(arr[i]), end="")
	
	if fmt == "float":
		print("{:.2f}/)".format(arr[-1]))
	elif fmt == "sci":
		print("{:.2e}/)".format(arr[-1]))

# Print out as fortran arrays
for s in ["W", "Fe"]:
	print(s)
	print_fortran_array(data[s]["Te"], "tes", "float")
	print_fortran_array(data[s]["ne"], "nes", "sci")

	# Now each 1-fprompt
	fnp = np.array([])
	for i in range(0, len(data[s]["ne"])):
		fnp = np.append(fnp, data[s]["fnp"][i])
	print_fortran_array(fnp, "fnp", "sci")
