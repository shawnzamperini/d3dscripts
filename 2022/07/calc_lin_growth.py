# Script to calculate the interchange turbelence linear growth rate.
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


"""
root = "/Users/zamperini/Documents/d3d_work/files/ts_"
shots = {190484:[3000, 4000]}
tss = []
for shot in shots.keys():

    # Load TS dictionary.
    with open("{}{}.pickle".format(root, shot), "rb") as f:
        ts = pickle.load(f)
        tss.append(ts)

    times = shots[shot]
    for time in times:

        # Find index closest to time.
        ts_times = ts["core"]["time"]
        tidx = np.argmin(np.abs(ts_times-time))
        print("{}: {}".format(time, ts_times[tidx]))

        # Pull out profiles during this time.
        ts["core"][]


    # Smooth the TS data some?
"""

# Load all the OMFIT profile data. The times within are only at the times
# of the RCP plunges.
root = "/Users/zamperini/Documents/d3d_work/files/OMFITprofiles_"
shots = [190440, 190442, 190484, 190485, 190486]
ncs = [netCDF4.Dataset("{}{}_FIT.nc".format(root, shot)) for shot in shots]
ncs_der = [netCDF4.Dataset("{}{}_DERIVED.nc".format(root, shot)) for shot in shots]
minor_radius = 0.584
rho_window = [0.95, 1.05]
mi = 931.49e6 / (3e8**2)
R = 1.94  # Core TS at same R

growth_rates = {}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))
for i in range(0, len(shots)):

    # Pull out the times and associated data.
    nc = ncs[i]
    shot = shots[i]
    times = nc["time"][:].data
    tes = nc["T_e"][:].data
    nes = nc["n_e"][:].data
    dnes = nc["dn_e_drho"][:].data
    rhos = nc["rho"][:].data

    growth_rates[shot] = []
    for j in range(0, len(times)):

        label = "{}.{}".format(shot, times[j])
        print()
        print(label)

        # Calculate radial gradient of ln(ne).
        ln_ne = np.log(nes[j])
        dln_ne = np.gradient(ln_ne, rhos*minor_radius)

        # Find the maximum (minimum since decreasing) dn/drho in the indicated
        # rho window.
        mask = np.logical_and(rhos>=rho_window[0], rhos<=rho_window[1])
        min_dne_idx = dln_ne[mask].argmin()
        #output[shot].append(dnes[j][min_dne_idx])
        min_dne_rho = rhos[mask][min_dne_idx]
        Ln = -dln_ne[mask][min_dne_idx]
        print("min_rho = {}".format(min_dne_rho))
        print("Ln = {}".format(Ln))

        # Calculate the sound speed at this location.
        fte = interp1d(rhos, tes[j])
        te = fte(min_dne_rho)
        cs = np.sqrt(2 * te / mi)
        print("te = {}".format(te))
        print("cs = {}".format(cs))

        # Calculate the interchange growth rate.
        growth = cs / np.sqrt(R * Ln)
        growth_rates[shot].append(growth)

        # Plots.

        ax1.plot(rhos, nes[j], label=label)
        ax2.plot(rhos, dln_ne, label=label)

#ax1.set_yscale("log")
ax1.grid(which="both")
ax2.legend()
fig.tight_layout()
fig.show()


lambdas = np.array([0.00568,
0.17578,
0.02618,
0.08170,
0.10178,
0.01092,
0.05127])
growth_rate = np.array([9.31E+03,
1.15E+04,
9.02E+03,
9.26E+03,
5.95E+03,
1.83E+04,
1.36E+04])
zeff = np.array([1.087,
1.066,
1.062,
1.182,
1.061,
1.154,
1.096])


from sklearn.linear_model import LinearRegression

X = np.array([growth_rate, zeff]).T
y = lambdas

model = LinearRegression().fit(X, y)

pred = model.predict(X)

fig, ax = plt.subplots()
ax.plot([0, 1], [0, 1], color="k")
ax.scatter(pred, y)
ax.grid()
fig.tight_layout()
fig.show()
