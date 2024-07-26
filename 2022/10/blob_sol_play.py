# Play around script to launch blobs and track them through the SOL, with the
# goal of getting average ne, Te values.
import numpy as np
from tqdm import tqdm
from scipy.stats import gengamma, norm
import matplotlib.pyplot as plt
import random
import sys
from scipy.optimize import fsolve, curve_fit
import matplotlib as mpl


# Inputs.
vr_mean    = 20
vr_std     = 375
vr_skew    = 0
gamma_a    = 1.8
gamma_c    = 0.80
nparts     = 20000
xpt_loc    = 15
qtim       = 1e-5
blob_ne    = 1e19
blob_te    = 100
blob_f     = 1e3  # Unused
seed_te    = 10  # Starting target Te.
niter      = 2
tau_te     = 0.0002
tau_ne     = 0.0005
blob_width = 0.01

# Constants.
mi = 931.49e6

# Here load in a MAFOT run that has the total field line lengths, and use that
# to create boundaries on our simulation grid.
# To-do.


# Simple 2D rectangular simulation volume. X = radial, Y = parallel.
maxx = 0.1
maxy = 100
xs = np.linspace(0, maxx, 100)
ys = np.linspace(0, maxy, 100)
X, Y = np.meshgrid(xs, ys)

# Generate a number of blobs to launch with vr's from the distribution.
launch_vrs = gengamma.rvs(gamma_a, gamma_c, loc=vr_mean, scale=vr_std, size=nparts)
below_zero = np.where(launch_vrs <= 0)
launch_vrs = np.delete(launch_vrs, below_zero)
print("Actual particles launched: {}".format(len(launch_vrs)))
mean_vr = launch_vrs.mean()

# Arrays for statistics.
counts = np.zeros(X.shape)
events = 0
mid_nes = np.zeros((niter, len(xs)))
mid_tes = np.zeros((niter, len(xs)))
vys = np.zeros(X.shape)
tetarg_warn = False

def exp_weight(t, tau):
    return np.exp(-t/tau)

# Follow blobs until deposition on the vessel.
for j in range(0, niter):
    print("Iteration: {}".format(j))

    # Calculate the parallel flow along each flux tube.
    prev_vys = vys.copy()
    vys = np.zeros(X.shape)
    for k in range(0, len(xs)):
        #if j == 0:
        #    prev_te = np.nan
        #else:
        #    prev_te = te[0,k]

        if j == 0:
            local_te = seed_te
        else:
            local_te = te[0,k]
            if np.isnan(local_te) or local_te == 0:
                #if np.isnan(prev_te) or prev_te == 0:

                # At this point the current iteration did not provide a new
                # Te, and the previous didn't either, so just use the
                # nearest Te value available. This will mainly apply in the
                # near-SOL where blobs don't generally reach. Elegant solution
                # from https://stackoverflow.com/questions/9537543/replace-nans-in-numpy-array-with-closest-non-nan-value
                mask = np.isnan(te[0])
                te[0, mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), te[0][~mask])
                local_te = te[0,k]

                    #local_te = seed_te
                #else:
                #    local_te = prev_te
                tetarg_warn = True
        local_cs = np.sqrt(2 * local_te / mi) * 3e8

        # Treat M as linear, increasing to cs at each target. Using point-slope.
        mid_idx = int(len(ys) / 2)
        slope = (local_cs) / (ys.max() - ys[mid_idx])
        local_vy = slope * (ys-ys.max()) + local_cs

        # Don't reassign if it's all zeros.
        if local_vy.sum() == 0:
            vys[:,k] = prev_vys[:,k]
        else:
            vys[:,k] = local_vy

    # Can we do a solver for the 1D fluid equations where ne, Te are now
    # prescribed? With nothing on the RHS:
    # d/ds(mi*n*v^2 + p||i + pe) = 0
    #tube_te = te[:,k]
    #tube_ne = ne[:,k]
    #tube_p = 2 * tube_te * tube_ne

    ne_weights = np.zeros(X.shape)
    te_weights = np.zeros(X.shape)
    ne_counts = np.zeros(X.shape)
    te_counts = np.zeros(X.shape)

    for i in tqdm(range(0, nparts)):

        # Choose starting y location between X-points.
        #y = float(xpt_loc + random.random() * (maxy - 2*xpt_loc))
        y = 0
        #while y < xpt_loc or y > (maxy-xpt_loc):
        while y <= 0 or y >= maxy:
            y = float(np.random.normal(maxy/2, maxy/5))

        x = 0.0
        vr = launch_vrs[i]
        t = 0.0
        dist = np.sqrt(np.square(X-x) + np.square(Y-y))
        loc_idx = np.where(dist==dist.min())

        # Record starting point in statistics.
        ne_weights[loc_idx] += 1.0
        ne_counts[loc_idx] += 1
        te_weights[loc_idx] += 1.0
        te_counts[loc_idx] += 1

        # Blob length is just the distance between the center and the nearest
        # X-point (and then an equal distance in the opposite direction). Later,
        # a gaussian will be "fit" into this length for weighing purposes.
        if y <= maxy / 2:
            blob_length = y - xpt_loc
        else:
            blob_length = (maxy - xpt_loc) - y

        while True:

            # Step sizes.
            vy = vys[loc_idx]
            ystep = vy * qtim
            xstep = vr * qtim
            #print(ystep)

            # Update location.
            y += ystep
            x += xstep
            t += qtim

            # Check for boundaries.
            if x >= maxx or y >= maxy:
                break

            # Update weights. No change for now.
            #ne_weight = exp_weight(t, tau_ne)
            #te_weight = exp_weight(t, tau_te)
            #ne_weight = 1.0
            #te_weight = 1.0
            ne_mult = exp_weight(t, tau_ne)
            te_mult = exp_weight(t, tau_te)

            # Tally its position.
            dist = np.sqrt(np.square(X-x) + np.square(Y-y))
            loc_idx = np.where(dist==dist.min())
            #ne_weights[loc_idx] += ne_weight
            #te_weights[loc_idx] += te_weight
            #ne_counts[loc_idx] += 1
            #te_counts[loc_idx] += 1
            events += 1

            # Blobs would have some sort of length and width. We can approximate
            # them as gaussian distributions in each direction, and then apply
            # the appropriately weighted values to the neighboring cells down to
            # specified minimum (to save computational time).
            # First find which radial bins the blob spans. Do this by seeing
            # how wide the blob is down to, e.g., 10% from it's center location.
            rad_coords = x - X[loc_idx[0][0]]
            rad_weights = norm.pdf(rad_coords, scale=blob_width)
            rad_weights = rad_weights / rad_weights.max()
            rad_cells = np.where(np.logical_and(rad_coords>=-blob_width/2, rad_coords<=blob_width/2))
            rad_weights = rad_weights[rad_cells]
            ne_weights[loc_idx[0][0]][rad_cells] += rad_weights * ne_mult
            ne_counts[loc_idx[0][0]][rad_cells] += 1
            te_weights[loc_idx[0][0]][rad_cells] += rad_weights * te_mult
            te_counts[loc_idx[0][0]][rad_cells] += 1

            # Next for the parallel, we will just assume they have constant
            # lengths which was chosen above. Only add weights and counts to
            # cells which the blob overlaps with.
            par_coords = y - Y[:,loc_idx[1]]
            par_cells = np.where(np.logical_and(par_coords>=-blob_length/2, par_coords<=blob_length/2))[0]
            par_weights = np.full(len(par_cells), 1.0).reshape(-1,1)
            ne_weights[:,loc_idx[1]][par_cells] += par_weights * ne_mult
            ne_counts[:,loc_idx[1]][par_cells] += 1
            te_weights[:,loc_idx[1]][par_cells] += par_weights * te_mult
            te_counts[:,loc_idx[1]][par_cells] += 1


    # Normalize counts and calculate plasma quantities. cell_volume cancels out,
    # but it doesn't necesarilly need to if the cell sizes were not all the same.
    #cell_volume = (xs[1] - xs[0]) * (ys[1] - ys[0]) * 1.0
    #nelec = cell_volume * blob_ne * blob_f * qtim
    #ne = ne_weights * nelec / nparts / cell_volume
    #ne = counts * blob_ne

    # Te is just the average Te of blobs counted within the cell. te_weights
    # contains the sum of the weights, and te_counts the number of blobs. E.g.
    # te_weights[i,j] = 2.5 and te_counts[i,j] = 3, then the average is 2.5 / 3.0 * blob_te.
    te = te_weights / te_counts * blob_te
    ne = ne_weights / ne_counts * blob_ne

    # Midplane profiles of ne and Te.
    mid_idx = np.argmin(np.abs(ys-ys.max()/2)) - 1
    r = X[mid_idx]
    mid_nes[j] = ne[mid_idx]
    mid_tes[j] = te[mid_idx]


# Target profiles.
targx = X[0]
targne = ne[0]
targte = te[0]

# Flux tube 1 cm from the separatrix.
tube_idx = np.argmin(np.abs(xs-0.03))
tube_s = Y[:,tube_idx]
tube_ne = ne[:,tube_idx]
tube_te = te[:,tube_idx]

# Grid of plots.
fig, axs = plt.subplots(3, 3, figsize=(9,8))
axs = axs.flatten()

#cmesh = axs[0].pcolormesh(X, Y, counts, shading="auto")
#cbar = fig.colorbar(cmesh, ax=axs[0])
#axs[0].set_xlabel("Radial (m)")
#axs[0].set_ylabel("Parallel (m)")

cmap = mpl.cm.get_cmap("plasma", 10)
nemax = np.nanmax(ne[1:-1, 1:-1])
cmesh = axs[0].pcolormesh(X, Y, ne, shading="auto", vmax=nemax, cmap=cmap)
cbar = fig.colorbar(cmesh, ax=axs[0])
axs[0].set_xlabel("Radial (m)")
axs[0].set_ylabel("Parallel (m)")
axs[0].axhline(xpt_loc, color="k")
axs[0].axhline(maxy-xpt_loc, color="k")

cmap = mpl.cm.get_cmap("inferno", 10)
temax = np.nanmax(te[1:-1, 1:-1])
cmesh = axs[1].pcolormesh(X, Y, te, shading="auto", vmax=temax, cmap=cmap)
cbar = fig.colorbar(cmesh, ax=axs[1])
axs[1].set_xlabel("Radial (m)")
axs[1].set_ylabel("Parallel (m)")
axs[1].axhline(xpt_loc, color="k")
axs[1].axhline(maxy-xpt_loc, color="k")

cmap = mpl.cm.get_cmap("coolwarm", 11)
cmesh = axs[2].pcolormesh(X, Y, vys, shading="auto", cmap=cmap)
cbar = fig.colorbar(cmesh, ax=axs[2])
axs[2].set_xlabel("Radial (m)")
axs[2].set_ylabel("Parallel (m)")

for i in range(0, mid_nes.shape[0]):
    axs[3].plot(r, mid_nes[i])
    axs[4].plot(r, mid_tes[i])
axs[3].set_xlabel("R (m)")
axs[3].set_ylabel("ne (m-3)")
axs[4].set_xlabel("R (m)")
axs[4].set_ylabel("Te (eV)")

axs[5].plot(targx, targne)
axs[5].set_xlim([0, maxx])
axs[5].set_xlabel("Distance from strike point (m)")
axs[5].set_ylabel("ne (m-3)")

axs[6].plot(targx, targte)
axs[6].set_xlim([0, maxx])
axs[6].set_xlabel("Distance from strike point (m)")
axs[6].set_ylabel("Te (eV)")

axs[7].plot(tube_s, tube_ne)
axs[7].set_xlabel("S (m)")
axs[7].set_ylabel("ne (m-3)")

axs[8].plot(tube_s, tube_te)
axs[8].set_xlabel("S (m)")
axs[8].set_ylabel("Te (eV)")

fig.tight_layout()
fig.show()
