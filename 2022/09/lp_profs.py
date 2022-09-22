import numpy as np
import matplotlib.pyplot as plt


# Generic inputs.
blob_peak = 50.0   # m from the target
blob_std = blob_peak / 2  # keep it mostly above X-point
xp_loc = 5.0  # in m from the target.
vr = 300  # m/s
A = 1  # The coefficient out front.
f = 3.0  # Flux expansion at the divertor.
nparts = 5000
te = 100.0

# Calculate sound speed.
mi = 931.49e6
cs = np.sqrt((2*te)/mi) * 3e8

# Choose initial starting point (x coordinate) from Guassian. Remove anything
# that is below the X-point.
launch_locs = np.random.normal(blob_peak, blob_std, nparts)
below_xp = np.where(launch_locs<= xp_loc)
launch_locs = np.delete(launch_locs, below_xp)
print("Actual particles launched: {}".format(len(launch_locs)))

targ_locs = np.zeros(len(launch_locs))
for i in range(0, len(launch_locs)):
    x = launch_locs[i]

    # Approximate length of blob as the distance between the launch location
    # and the X-point.
    #blob_length = x - xp_loc

    # Define point at which it strike the target as point where the tip of the
    # blob hits.
    y = vr / cs * x * f
    targ_locs[i] = y

# Bin data into histogram.
bin_locs = np.histogram(targ_locs, 100)
plotx = np.array([(bin_locs[1][i]+bin_locs[1][i+1])/2 for i in range(0, len(bin_locs[1])-1)])

fig, ax1 = plt.subplots()

ax1.plot(plotx*100, bin_locs[0])
ax1.set_xlabel("Distance along target (cm)", fontsize=14)
ax1.set_ylabel("ne (m-3)", fontsize=14)

fig.tight_layout()
fig.show()
