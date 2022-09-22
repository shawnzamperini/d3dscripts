import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle
from scipy.stats import skewnorm

# As a test case, let's try and match a particularly nice LP profile.
with open("lp174175.pickle", "rb") as f:
    data = pickle.load(f)

# Generic inputs.
blob_peak = 50.0   # m from the target
blob_std = blob_peak / 2  # keep it mostly above X-point
xp_loc = 5.0  # in m from the target.
vr_mean = 200  # m/s
vr_std = 200
vr_skew = 0.5
A = 3e17  # The coefficient out front.
f = 3.0  # Flux expansion at the divertor.
nparts = 5000
tesep = 50.0
lambda_te = 0.1
rs = np.linspace(0, 0.1, 100)
qtim = 1.0e-7

# Calculate sound speed.
mi = 931.49e6
te = tesep * np.exp(-rs/lambda_te)
cs = np.sqrt((2*te)/mi) * 3e8

# Choose initial starting point (x coordinate) from Guassian. Remove anything
# that is below the X-point.
#launch_locs = np.random.normal(blob_peak, blob_std, nparts)
#below_xp = np.where(launch_locs <= xp_loc)
#launch_locs = np.delete(launch_locs, below_xp)
#print("Actual particles launched: {}".format(len(launch_locs)))

# Choose blobs vr's. Anything below zero remove.
#launch_vrs = np.random.normal(vr_mean, vr_std, nparts)
launch_vrs = skewnorm.rvs(a=vr_skew, loc=vr_mean, scale=vr_std, size=nparts)
below_zero = np.where(launch_vrs <= 0)
launch_vrs = np.delete(launch_vrs, below_zero)
print("Actual particles launched: {}".format(len(launch_vrs)))

targ_locs = np.zeros(len(launch_vrs))
for i in tqdm(range(0, len(launch_vrs))):

    # Get starting launch location.
    #x = launch_locs[i]

    # Assume all blobs extend to X-point, so that's the distance they need to
    # travel before hitting the target.
    x = xp_loc


    # Follow particle in 1D sense.
    xpos = x
    rpos = 0.0
    vr = launch_vrs[i]
    while True:

        # Perform radial step, see what nearest cs is (which we treat as v parallel).
        rpos = rpos + vr * qtim
        ridx = np.argmin(np.abs(rs-rpos))
        vpar = cs[ridx]

        # Perform parallel step, negative because towards target.
        xpos = xpos - vpar * qtim

        # If xpos goes below zero, we've hit the target. Record location. The
        # flux expansion factor, f, is dubious in this simple approximation.
        if xpos < 0:
            targ_locs[i] = rpos * f
            break

    # Approximate length of blob as the distance between the launch location
    # and the X-point.
    #blob_length = x - xp_loc

    # Define point at which it strike the target as point where the tip of the
    # blob hits.
    #y = vr / cs * x * f
    #targ_locs[i] = y

# Bin data into histogram.
bin_locs = np.histogram(targ_locs, 100)
plotx = np.array([(bin_locs[1][i]+bin_locs[1][i+1])/2 for i in range(0, len(bin_locs[1])-1)])

# Also bin the vr-pdf.
plot_vr = np.linspace(0, vr_mean*3)
pdf = skewnorm.pdf(plot_vr, a=vr_skew, loc=vr_mean, scale=vr_std)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

ax1.scatter(data["rmrs"]*100, data["ne"], c=data["color"], zorder=5)
ax1.plot(data["rmrs"]*100, data["nefilt"], color="k", zorder=10)
ax1.set_xlabel("R-Rsep (m)", fontsize=14)
ax1.set_ylabel("ne (m-3)", fontsize=14)

ax1.plot(plotx*100, A*bin_locs[0], color="r", lw=2, zorder=20)
#ax1.set_ylabel("Counts", fontsize=14)

ax2.plot(plot_vr, pdf)
ax2.set_xlabel("vr (m/s)", fontsize=14)
ax2.set_ylabel("PDF", fontsize=14)

fig.tight_layout()
fig.show()
