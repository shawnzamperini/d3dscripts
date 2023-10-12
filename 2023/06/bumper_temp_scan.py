import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from tqdm import tqdm

# Some inputs.
shot_start = np.linspace(30, 34, 75)  # Maximum this can be is 35 (for a 5 second shot). Minimum anything less than 30 doesn't make sense.
shot_duration = np.linspace(4.5, 5.5, 5)  # It's pretty much 5 seconds but can add some play.
# shot_duration = [5]
heat_flux = np.linspace(1e6, 10e6, 5)  # A pretty intense but reasonable range to scan in.

# The following slopes and intercepts are actually for 1/value.
therm_diff_slope = np.linspace(25, 45, 5)  # CFC is 37.
# therm_cond_slope = np.linspace(5e-6, 2e-5, 5)  # Values from source are 6.56e-6 and 1.007, 1.16 and 1.19e-5.
therm_cond_slope = [4.48e-6]  # We have this data for the CFC of the limiters.
# therm_cond_intc = np.linspace(0.005, 0.01, 5)  # Really just eyeballing this from the plot.
therm_cond_intc = [1.02E-2]  # We have this data for the CFC of the limiters.


# Assemble all possible combinations of values.
# combos = np.stack(np.meshgrid(shot_start, shot_duration, therm_diff_slope, therm_cond_slope, therm_cond_intc, heat_flux), -1).reshape(-1, 6)

# Alternatively just assemble a random combination of a set number of combinations in the range of the values.
def assemble_rans(arr, nvals):
    if len(arr) == 1:
        return np.full(nvals, arr[0])
    else:
        return arr.min() + (arr.max() - arr.min()) * np.random.random(nvals)


combos = np.stack(
    np.meshgrid(assemble_rans(shot_start, len(shot_start)),
                assemble_rans(shot_duration, len(shot_duration)),
                assemble_rans(therm_diff_slope, len(therm_diff_slope)),
                assemble_rans(therm_cond_slope, len(therm_cond_slope)),
                assemble_rans(therm_cond_intc, len(therm_cond_intc)),
                assemble_rans(heat_flux, len(heat_flux))),
    -1).reshape(-1, 6)

# Just hardcoding in a sample peak of high temperature. Convert temperature from C to K.
lim_temp = np.array(
    [254.9500122, 253.607193, 252.2949066, 255.7129822, 648.5786133, 612.017334, 554.5508423, 518.0201416,
     490.7975769, 472.9442139, 458.9361877, 449.4448853, 441.5100708])
lim_temp = lim_temp + 273.0
lim_time = np.arange(0, 10 * len(lim_temp), 10)


def therm_diff(T, slope=52):
    """
    Thermal diffusivity in [m2/s]. Defaults for "Analysis 6510".
    """
    # print("cond: {:}".format(1 / slope * T))
    return 1 / (slope * T)


def therm_cond(T, slope=1.007e-5, intc=6.8008e-3):
    """
    Thermal conductivity in [W/K*m]. Defaults for "Analysis 6510".
    """
    return 1 / (slope * T + intc)


def surf_temp(temp, t, qw, tp, temp0, tds, tcs, tci):
    """
    Calculate the surface temperature of graphite. If t < tp we are still in the heating phase and use the correpsonding
    equation, otherwise use the cooling down equation.

    temp (float): Surface temperature, what we are solving for.
    t (float): Time in s.
    qw (float): Impinging heat flux in W/m2.
    tp (float): The pulse length at which the heat flux turns off in s.
    temp0 (float): The starting temperature of the surface.
    """

    if t < tp:

        # Rewritten to eliminate c in favor of alpha (and the density cancels out).
        return 2 * qw * np.sqrt(t) / (
            np.sqrt(np.pi * therm_cond(temp, tcs, tci) ** 2 / therm_diff(temp, tds))) + temp0 - temp
        # return qw * np.sqrt(t) / 10000
    else:
        return 2 * qw * np.sqrt(therm_diff(temp, tds)) / (therm_cond(temp, tcs, tci) * np.sqrt(np.pi)) * (np.sqrt(t)
                                                                                                          - np.sqrt(
                    t - tp)) + temp0 - temp
        # return qw / 10000 * (np.sqrt(t) - np.sqrt(t - tp)) + temp0 - temp


# # Test.
# qw = 7.0e6
# tp = 5.0
# temp0 = lim_temp[3]
# tend = 40.0
# nvals = 200
# temps = np.zeros(nvals)
# times = np.linspace(0, tend, nvals)
# prev_temp = temp0
# #for i in range(0, nvals):
#     #root = fsolve(surf_temp, prev_temp, args=(times[i], qw, tp, temp0))
# #    res = least_squares(surf_temp, prev_temp, bounds=(1, np.inf), args=(times[i], qw, tp, temp0))
#     #print("{:5.2f}  {:7.2f}".format(times[i], res.x[0]))
# #    temps[i] = res.x[0]
# #    prev_temp = res.x[0]

# Scan through every possible value combination and save the ones where the end temperature is within a set percentage
# of what the corresponding final value is (i.e., lim_temp[4]).
count = 0
results = {}
nvals = 100
temp0 = lim_temp[3]
temp_end = lim_temp[4]
temp_end_t = lim_time[4]
tend = 20
times = np.linspace(0, tend, nvals)
margin = 0.10

# Parallelize this???
k = 0
for combo in tqdm(combos):

    # Extract all the variables.
    t0 = combo[0]
    tp = combo[1]
    tds = combo[2]
    tcs = combo[3]
    tci = combo[4]
    qw = combo[5]

    # Setup arrays.
    temps = np.zeros(nvals)

    # Loop through for each time.
    prev_temp = temp0
    for i in range(0, nvals):
        res = least_squares(surf_temp, prev_temp, bounds=(1, np.inf), args=(times[i], qw, tp, temp0, tds, tcs, tci))
        temps[i] = res.x[0]
        prev_temp = res.x[0]

    # See what the final predicted temperature at the time of the first measurement after the pulse is.
    # need to use t0 here as well as temp_end_t to  find the nearest time in times to get the temperature for.
    combo_times = t0 + times  # The times offset to align with the limiter data.
    idx = np.abs(combo_times - temp_end_t).argmin()
    if np.abs(temps[idx] - temp_end) / temp_end <= margin:
        # Within margin of error, save the profile.
        # print("match! {}".format(i))
        results[k] = {"t0": t0, "tp": tp, "tds": tds, "tcs": tcs, "tci": tci, "qw": qw, "times": combo_times,
                      "temps": temps}

    k += 1

# Extract the average/median max temperatures.
tmp_temps = []
for i in results.keys():
    tmp_temps.append(results[i]["temps"].max())
avg_max = np.mean(tmp_temps)
std_max = np.std(tmp_temps)
med_max = np.median(tmp_temps)
print("Peak Limiter Temperature")
print(" Mean: {:} +/- {:} K".format(int(avg_max), int(std_max)))
print(" Median: {:} K".format(int(med_max)))

fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.plot([0], [0], color="tab:red", lw=2, label="1D Heat Flux Est.")
for i in results.keys():
    ax1.plot(results[i]["times"], results[i]["temps"], color="tab:red", lw=1, alpha=0.05)
ax1.scatter(lim_time, lim_temp, zorder=15, marker="^", edgecolors="k", color="tab:red", label="TC Data", s=75)
ax1.set_xlim([25, 55])
ax1.set_ylim([400, 2000])
ax1.grid(alpha=0.3, zorder=5)
ax1.set_xlabel("Time, arbitrary (s)", fontsize=12)
ax1.set_ylabel("Limiter temperature (K)", fontsize=12)
ax1.legend()
fig.tight_layout()
fig.show()
