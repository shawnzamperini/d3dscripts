import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

path = "/Users/zamperini/My Drive/Research/Documents/2023/06/bumper_tile_temps.xlsx"
df = pd.read_excel(path)

# Get subset of the highest temperature peak.
subset = df.iloc[575:588]

# Get X value in just seconds from the start of the data.
ts = (subset["t_stamp"] - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
time = ts - ts.min()

# Y data is just the temperature (in seconds).
temp = subset["Limiter310_4"]

# Igor says the temperature evolves as the square root of time for an infinite square pulse, which we treat the shot
# to be. Thus, to estimate the possible bumper temperature at the end of the shot (which we take to be the max), we fit
# a function of sqrt(t) to the data after the peak and back track to see what it could've been.
after_time = time[5:]
after_temp = temp[5:]

def sqrt_time(t, a, b):
    return a * np.sqrt(t) + b

offset = 15
popt, pcov = curve_fit(sqrt_time, after_time-offset, after_temp)
time_fit = np.linspace(50-offset, 90-offset, 100)
temp_fit = sqrt_time(time_fit, *popt)
time_fit2 = np.linspace(40-offset, 50-offset, 20)
temp_fit2 = sqrt_time(time_fit2, *popt)


# Plot.
fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.fill_between([40, 50], [100, 100], [900, 900], color="tab:purple", alpha=0.3)
ax1.plot(time_fit+offset, temp_fit, label=r"$\mathdefault{\sqrt{t}}$ Fit", lw=3, color="tab:red")
ax1.plot(time_fit2+offset, temp_fit2, lw=3, color="tab:red", linestyle="--")
ax1.plot(time, temp, marker="^", label="Data", lw=0, ms=7, mec="k", color="tab:red")
ax1.text(51, 725, "#195128", color="tab:purple", fontsize=12)
ax1.set_xlabel("Time, arbitrary (s)", fontsize=12)
ax1.set_ylabel("Temperature (C)", fontsize=12)
ax1.set_xlim(30, 110)
ax1.set_ylim(200, 800)
ax1.legend()
ax1.grid(alpha=0.3)
ax1.set_title("4/10/23 Limiter 310-4")
fig.tight_layout()
fig.show()
