# An Attempt to simulate reerosion from a deposition profile from 3DLIM.
import numpy as np
import random
import matplotlib.pyplot as plt
from tqdm import tqdm


# Inputs.
test_exp           = True
reerosion_strength = 1.0
launch_energy      = 5
mass_ion           = 12   # C = 12, Si = 28
charge_ion         = 1
te                 = 5
ti                 = 2*te
ne                 = 1e18
#vz_mult            = 0.9
rad_vel            = 375  # m/s
delta_t            = 1e-8
precision          = 2         # Precision of which to convert deposition units to.
n_parts            = None    # If None just do all of them.
ff_mod             = 1

# Mass of impurity in kg.
mass_ion_kg = mass_ion * 1.66E-27

# If test_exp then don't bother with real data and use an exponential.
if test_exp:
    dep_x = np.linspace(0, 0.1, 100)
    dep_y = np.exp(- dep_x / 0.01)
    reerosion_zone = [0, 0.02]
    #reerosion_cdf = np.cumsum(dep_y / dep_y.sum())
    #reerosion_zone = [dep_x.min(), dep_x.max()]

# Normalize the data such that you can go from floats --> ints. Do this by
# rounding to the nearest 6 digits and multiplying by 1e6 (e.g. any bin in X
# location counts particles down to millionth precision).
dep_y = np.array(np.round(dep_y, precision) * np.power(10, precision), dtype=int)
net_y = np.copy(dep_y)

# Particles within the reerosion zone get reeroded.
zone_mask = np.logical_and(dep_x >= reerosion_zone[0], dep_x <= reerosion_zone[1])
if type(n_parts) != type(None):
    n_reerode = n_parts
else:
    n_reerode = int(dep_y[zone_mask].sum() * reerosion_strength)

# The parallel force will just be a constant Ffric at a single te.
col_log = 15
tau_s = 1.47E13 * mass_ion * ti * np.sqrt(ti / 2.0) / \
        ((1 + 2.0 / mass_ion) * ne * np.power(charge_ion, 2) * col_log)
vi = -np.sqrt((te + ti) * 1.609e-19 / (2.0 * 1.66e-27))   # Just cs, negative bc it's "downwards" towards the probe surface.
#vz = vz_mult * vi
#ff = mass_ion_kg * (vi - vz) / tau_s
#accel = ff / mass_ion_kg

# Just launch at same velocity every time for now.
launch_v = np.sqrt(2 * launch_energy * 1.609e-19 / mass_ion_kg)

# Arrays for stat keeping.
x_dists = np.zeros(n_reerode)
erodes = np.zeros(len(dep_x))

# Track one reeroded particle at a time.
for i in tqdm(range(0, n_reerode)):

    # Location of reeroded particle chosen at random, find nearest location in
    # array. Count as reeroded (negative particle), don't go below zero.
    erode_x = random.uniform(reerosion_zone[0], reerosion_zone[1])
    erode_idx = np.argmin(np.abs(dep_x - erode_x))
    #erode_idx = np.argmin(np.abs(reerosion_cdf - random.random()))
    #erode_x = dep_x[erode_idx]
    if net_y[erode_idx] == 0:
        continue
    else:
        net_y[erode_idx] -= 1
    erodes[erode_idx] += 1

    # Let's just use DIVIMP launch option #3 because it's easy.
    # Vel/angle flag 0 : theta =+/-asin($), $ in (0,1)
    launch_ang = (-1 + random.random() * 2) * np.arcsin(random.random())

    # Distance along the probe is just a simple kinematic equation.
    vy_int = launch_v * np.cos(launch_ang)

    start_x = erode_x

    # The X velocity is chosen from a Gaussian.
    vx = random.gauss(rad_vel, 10)

    # Track particle until it redeposits (y = 0).
    x = start_x
    y = 0.001  # Start 1 mm off of surface I guess.
    vy = vy_int
    count = 0
    while True:

        # X and Y steps are calculated up front.
        x_step = vx * delta_t
        y_step = vy * delta_t

        # If y has gone below zero then it's redepositied.
        if y <= 0:
            redep_idx = np.argmin(np.abs(dep_x - x))
            break

        # Update values for next iteration.
        x += x_step
        y += y_step
        ff = mass_ion_kg * (vi - vy) / tau_s * ff_mod
        vy = vy + (ff / mass_ion_kg) * delta_t
        count += 1

    # Accounting for where it landed.
    net_y[redep_idx] += 1
    x_dists[i] = x - start_x

print("Average re-eroded distance: {:.2} cm".format(x_dists.mean() * 100))

plt.rcParams['font.family'] = 'Century Gothic'
fig, ax1 = plt.subplots()
ax1.plot(dep_x*100, dep_y, color="k", label="Inital", lw=2)
ax1.plot(dep_x*100, net_y, color="k", linestyle="--", label="With re-erosion", lw=2)
ax1.set_xlabel("Distance along probe (cm)", fontsize=14)
ax1.set_ylabel("Deposition (arbitrary units)", fontsize=14)
ax1.grid()
ax1.tick_params(labelsize=14)
ax1.legend(fontsize=14, facecolor="white")
fig.tight_layout()
fig.show()
