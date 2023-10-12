import matplotlib.pyplot as plt
import numpy as np


meas_s = 80
smax   = 100
meas_mach = 0.15

# s=0 is the inner target.
s = np.linspace(0, smax, 100)

# Select correct slope for which side of the pyramid we're on to calculate stagnation point.
slope2 = (meas_mach - 1) / (meas_s - smax)
stag_s = -1 / slope2 + smax
slope1 = 1 / stag_s

mach = np.zeros(len(s))
mach[s < stag_s] = slope1 * s[s < stag_s] - 1
mach[s > stag_s] = slope2 * (s[s > stag_s] - smax) + 1

fig, ax = plt.subplots(figsize=(4, 3))
ax.axhline(0, color="k")
ax.plot(s, mach, color="tab:red", zorder=10, lw=3)
ax.plot([0, smax], [-1, 1], color="tab:red", zorder=11, linestyle="--", lw=3)
ax.scatter(meas_s, meas_mach, marker="x", s=100, color="k", zorder=15, linewidths=2)
ax.grid(alpha=0.3)
ax.set_xlabel("s (m)", fontsize=12)
ax.set_ylabel("Mach Number", fontsize=12)
fig.tight_layout()
fig.show()