# This program is a Monte-Carlo simulation of blobs travelling through a slab
# geometry. The blobs are launched from the core portion of the separatrix
# from a normal radial velocity distribution. They keep the same radial
# velocity until they deposit on the boundaries. The parallel speed is set to
# the sound speed, which is calculated for each cell based off a basic Te
# prescription.
import numpy as np
from tqdm import tqdm
import pickle


# Inputs.
nblobs     = 300
te_sep     = 100
te_decay   = 0.05
fblob_sep  = 2e3
sim_width  = 0.15
sim_length = 100
pfz_length = [5, sim_length/2]   # Applied to both ends.
time_step  = 1e-7
vr_sep     = 2000  # Mean of distribution.
vr_width   = 500  # Width of distribution.
vr_wall    = 0.1  # Will linearly decrease vr down to this fraction at the wall.
cs_mult    = 10
vy_mode    = "constant"  # One of "cs" or "constant"
vy_const   = 300000


# First setup the grid. Simple rectangular mesh is all it is.
# X = radial, Y = parallel, coordinates are the center of the bin.
xbins = np.linspace(0, sim_width, 101)
ybins = np.linspace(0, sim_length, 101)
xcenter = [xbins[i] + (xbins[i+1] - xbins[i])/2 for i in range(0, len(xbins[:-1]))]
ycenter = [ybins[i] + (ybins[i+1] - ybins[i])/2 for i in range(0, len(ybins[:-1]))]
X, Y = np.meshgrid(xcenter, ycenter)

# Simplest SOL is a sheath-limited, which means just set Te to constant
# value. ne not actually needed for the basic sim.
te = te_sep * np.exp(-X / te_decay)

# Setup sound speed to linearly increase towards each target.
if vy_mode == "cs":
    md = 2 * 931.494 * 1e7 / 3e8**2
    cs_fact = 2 * (Y - 0) / (sim_length - 0) - 1
    cs = np.sqrt(2 * te / md) * cs_fact * cs_mult

# Just use a constant velocity with a stepwise change halfway between.
elif vy_mode == "constant":
    cs_fact = np.zeros(Y.shape)
    for i in range(0, Y.shape[0]):
        for j in range(0, Y.shape[1]):
            if Y[i, j] > sim_length / 2:
                cs_fact[i, j] = 1
            else:
                cs_fact[i, j] = -1
    cs = np.full(Y.shape, vy_const) * cs_fact

# Arrays to count deposition of blobs.
targ1 = X[0]
targ2 = X[0]
wall  = Y[:,0]
targ1_dep = np.zeros(targ1.shape)
targ2_dep = np.zeros(targ2.shape)
wall_dep  = np.zeros(wall.shape)

# Particle arrays.
nsteps = np.zeros(nblobs, dtype=int)
launch_ys = np.zeros(nblobs)
blob_dens = np.zeros(X.shape)

# Launch one particle at a time, following it until it hits one of the boundaries.
for i in tqdm(range(0, nblobs)):

    # Launch uniformly between the start and end of the core region at X = 0
    y = np.random.uniform(pfz_length[0], sim_length-pfz_length[1])
    launch_ys[i] = y

    # Normal distributed velocity.
    vx = 0
    while vx <= 0:
        vx = np.random.normal(vr_sep, vr_width)

    # Setup vr array where vr linearly decreases towards the wall (point-slope).
    m = (vx - vr_wall*vx) / (0 - sim_width)
    vxs = m * X[0] + vx

    # Track particle until it hits one of the bounds.
    t = 0; x = 0; end = False; nstep = 0
    while not end:

        # Find current cell, get parallel velocity.
        dist = np.sqrt(np.abs(X-x) + np.abs(Y-y))
        idx = np.where(dist == dist.min())
        vy = cs[idx]
        vx = vxs[idx[1][0]]

        # Perform step, update coordinates.
        x = x + vx * time_step
        y = y + vy * time_step
        nstep += 1

        # Check if outside boundaries. If so, score what target bin it landed in.
        out = np.logical_or(x < 0, x > sim_width)
        if out:
            dist = np.abs(y - wall)
            idx = np.where(dist == dist.min())
            wall_dep[idx] += 1

            nsteps[i] = nstep
            break
        out = np.logical_or(y < 0, y > sim_length)
        if out:

            if y < 0:
                dist = np.abs(x - targ1)
                idx = np.where(dist == dist.min())
                targ1_dep[idx] += 1
            else:
                dist = np.abs(x - targ1)
                idx = np.where(dist == dist.min())
                targ2_dep[idx] += 1

            nsteps[i] = nstep
            break

        # Score it in arrays.
        blob_dens[idx] += 1

# Pickle output for plotting routine.
output = {"launch_ys":launch_ys, "blob_dens":blob_dens, "X":X, "Y":Y, "te":te,
    "targ1_dep":targ1_dep, "targ2_dep":targ2_dep, "wall_dep":wall_dep,
    "fblob_sep":fblob_sep, "nblobs":nblobs, "targ1":targ1, "targ2":targ2}

with open("slab_blob_output.pickle", "wb") as f:
    pickle.dump(output, f)
