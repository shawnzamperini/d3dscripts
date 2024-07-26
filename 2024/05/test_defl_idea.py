import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


if __name__ == "__main__":
    show_output = True
else:
    show_output = False

# Constants for W5+ and D+
amu = 1.660e-27  # kg
ev = 1.602e-19  # C
eps0 = 8.85e-12
mi = amu * 2.014  # kg
qi = ev  # C
me = 9.10938188e-31  # kg
#mz = amu * 183.38 # kg
mz = me
qz = -ev  # C
qe = -ev  # C

# Some simulation parameters
te = 1 * ev  # J
ne = 1e17 # m-3
ni = ne
#vz = 1e3  # m/s
vz = np.sqrt(2 * te / mz)
ignore_e_coll = False

# Reduced masses
mr_zi = mi * mz / (mi + mz)
mr_ze = me * mz / (me + mz)

# Number of pairs of D+ and e the impurity encounters over dt
lambda_d = np.sqrt(eps0 * te / (ne * np.square(ev))) 
n = int(np.pi * ni * np.power(lambda_d, 3))
#n = 100
#n = int(np.pi * eps0 * te * dt * vz / (ev * ev))
dt = lambda_d / vz
if show_output:
    print("Number of pairs: {}".format(n))
    #dt = 1e-8  # s
    print("Effective dt: {:.2e}".format(dt))

# Grab ion velocities from a Maxwellian distribution.
ti = te
vi_avg = np.sqrt(2 * ti /  mi)
vis = np.random.normal(vi_avg, vi_avg / 4, size=n) 
#vis = np.full(n, vi_avg)
vis[vis < 0] = 0.0
ve_avg = np.sqrt(2 * te / me)
ves = np.random.normal(ve_avg, ve_avg / 4, size=n)
ves[ves < 0] = 0.0

# Plot to make sure velocities are reasonable. 
if show_output:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    ax1.hist(vis, bins=15)
    ax2.hist(ves, bins=15)
    ax1.set_xlabel("vi (m/s)")
    ax2.set_xlabel("ve (m/s)")
    ax1.set_ylabel("Counts")
    fig.tight_layout()
    fig.show()

# Calculate the CoM initial velocity
v0s_i = np.sqrt(np.square(vz) + np.square(vis))
v0s_e = np.sqrt(np.square(vz) + np.square(ves))

# The impact parameters will be randomly distributed between 0 and the 
# Debye length.
#bs_i = np.random.random(n) * lambda_d
#bs_i = (2 * np.random.random(n) - 1) * lambda_d  # Equal probability of being "above" or "below" ion for collision.
#bs_e = np.random.random(n) * lambda_d

# Generate impacts parameters from a 1/r distribution. Need to define a non-zero minimum value.
min_b = 1e-8
max_b = lambda_d
bs_i = max_b - min_b * np.exp(np.random.random(n) * np.log(max_b / min_b))
bs_e = max_b - min_b * np.exp(np.random.random(n) * np.log(max_b / min_b))
#bs_i = np.full(n, lambda_d)
#fig, ax = plt.subplots(figsize=(5, 4))
#hist = ax.hist(bs_i, bins=25)
#ax.set_xlabel("Impact parameter (m)")
#fig.tight_layout()
#fig.show()

# Calculate a constant term ahead of time so we don't waste computational time.
#defl_angs_i = 2 * np.arctan2(1, (4 * np.pi * eps0 * mr_zi / (qi * qz) * np.square(v0s_i) * bs_i))
defl_angs_i = 2 * np.arctan(1 / (4 * np.pi * eps0 * mr_zi / (qi * qz) * np.square(v0s_i) * bs_i))
if ignore_e_coll:
    defl_angs_e = np.zeros(n)
else:
    defl_angs_e = 2 * np.arctan(1, (4 * np.pi * eps0 * mr_ze / (qe * qz) * np.square(v0s_e) * bs_e))

#defl_angs_i = np.array([v if v < np.pi/2 else 2*np.pi - v for v in defl_angs_i])

# Calculate all the deflection calculations then sum up.
tot_defl_ang_i = 0
tot_defl_ang_e = 0
#defl_angs_i = []
#defl_angs_e = []
defl_freq_90_count = 0  # Increment everytime we deflect 90 degrees.
defl_freq_90_tracker = 0.0
for i in tqdm(range(n)):

    # Deflection from deuterium
    #defl_i = 4 * np.pi * eps0 * mr_zi / (qi * qz) * np.square(v0s_i[i]) * bs_i[i]
    #defl_ang_i = 2 * np.arctan(1 / defl_i)
    #tot_defl_ang_i += defl_ang_i

    # Deflection from electrons. Option to ignore so we can make comparisons to 
    # theory. 
    #if ignore_e_coll:
    #    defl_e = 0.0
    #    defl_ang_e = 0.0
    #else:
    #    defl_e = 4 * np.pi * eps0 * mr_ze / (qe * qz) * np.square(v0s_e[i]) * bs_e[i]
    #    defl_ang_e = 2 * np.arctan(1 / defl_e)
    #tot_defl_ang_e += defl_ang_e

    # Append to list for plotting later.
    #defl_angs_i.append(defl_ang_i)
    #defl_angs_e.append(defl_ang_e)

    # Keep track of everytime we deflect 90 degrees so we can calculate a frequency later.
    defl_freq_90_tracker += (defl_angs_i[i] + defl_angs_e[i])
    if np.abs(defl_freq_90_tracker) > (np.pi / 2):
        defl_freq_90_tracker = 0
        defl_freq_90_count += 1

#defl_angs_i = np.array(defl_angs_i)
#defl_angs_e = -np.array(defl_angs_e)  # Minus for plotting, just remember this?

# Downsize to 0-2pi.
#tot_defl_ang_i = tot_defl_ang_i % (2 * np.pi)
#tot_defl_ang_e = tot_defl_ang_e % (2 * np.pi)
tot_defl_ang_i = np.sum(defl_angs_i) % (2 * np.pi)
tot_defl_ang_e = np.sum(defl_angs_e) % (2 * np.pi)
tot_defl_ang = (tot_defl_ang_i + tot_defl_ang_e) % (2 * np.pi)

if show_output:
    print("Total deflection: {:.2f} rad ({:.2f} deg)".format(tot_defl_ang, np.degrees(tot_defl_ang)))
    print("  Deuterium: {:.2f} radians   Average: {:.2e}".format(tot_defl_ang_i, np.mean(defl_angs_i)))
    print("  Electron:  {:.2f} radians   Average: {:.2e}".format(tot_defl_ang_e, np.mean(defl_angs_e)))
    print("*******************************")
    print("Total Average: {:.2e}".format(np.mean(defl_angs_i) + np.mean(defl_angs_e)))
    print("*******************************")

# See how often we hit 90 degrees. 
#defl_90_freq = defl_freq_90_count / dt
#defl_90_freq_i = defl_angs_i.sum() / dt * np.pi / 2.0
#defl_90_freq_e = defl_angs_e.sum() / dt * np.pi / 2.0
defl_90_freq = (defl_angs_i.sum() + defl_angs_e.sum()) / dt * np.pi / 2.0
if show_output:
    print("90 degree deflection frequency: {:.2e} Hz  ({} counts)".format(defl_90_freq, defl_freq_90_count))

# Considering this average approach...
avg_freq_90_i = (np.mean(defl_angs_i)/dt) / (np.pi/2)
std_freq_90_i = (np.std(defl_angs_i)/dt) / (np.pi/2)
if show_output:
    print("----------")
    print("Average 90 degree deflection frequency: {:.2e} +/- {:.2e}".format(avg_freq_90_i, std_freq_90_i))

# Calculate the value from Freidberg for comparison. 
ni = ne
lnalpha_fact = np.log(np.power(eps0, 3/2) * 4 * np.pi * amu / np.power(ev, 5/2))
vtz = np.sqrt(2 * ti / mz)
lnalpha_corr2 = lnalpha_fact + np.log(np.sqrt(te / ev / ne) * np.square(vtz) * (mr_zi / amu) / np.abs(qz/ev * qi/ev))
nu_fr_fact = np.power(ev, 4) * np.power(amu, 3/2) / (4 * np.pi * np.square(eps0) * np.square(amu) * np.power(2 * ev, 3/2))
nu_fr_corr2 = nu_fr_fact * np.square(qz/ev) * np.square(qi/ev) * ni * lnalpha_corr2 / ((mz/amu) * (mr_zi/amu)) / (np.power((ti/ev) / (mz/amu), 3/2) + 1.3 * np.power((ti/ev) / (mi/amu), 3/2))
if show_output:
    print("Freidberg Z-i 90 collision frequency: {:.2e} (every {:.2e} s)".format(nu_fr_corr2, 1/nu_fr_corr2))


# Plot to see what average deflections are.
if show_output:
    plt_defl_angs_i = np.abs(defl_angs_i)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    for ax in (ax1, ax2):
        for rad in np.arange(0, 2 * np.pi, np.pi / 2):
            ax.axvline(rad, color="k", linestyle="--")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim([1e-6, 2 * np.pi])
        ax.set_xlabel("|Deflection angle|")
    hist1 = ax1.hist(plt_defl_angs_i, bins=np.geomspace(np.min(plt_defl_angs_i), np.max(plt_defl_angs_i), 50), color="tab:red")
    if not ignore_e_coll:
        hist2 = ax2.hist(defl_angs_e, bins=np.geomspace(np.min(defl_angs_e), np.max(defl_angs_e), 25), color="tab:red")
    ax1.axvline(np.mean(plt_defl_angs_i), color="k", lw=5)
    ax1.axvline(np.mean(plt_defl_angs_i), color="tab:red", lw=3)
    ax2.axvline(np.mean(defl_angs_e), color="k", lw=5)
    ax2.axvline(np.mean(defl_angs_e), color="tab:red", lw=3)
    ax1.set_title("Deuterium")
    ax2.set_title("Electrons")
    fig.tight_layout()
    fig.show()

