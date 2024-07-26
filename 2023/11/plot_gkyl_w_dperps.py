import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Entries are for W1+, 5+, 10+, 15+.
xlpath = "/Users/zamperini/My Drive/Research/Documents/2023/11/gkyl_imp_dperps.xlsx"
xl = pd.read_excel(xlpath, skiprows=1)

# W1+ and W15+ are the bounding curves.
r1 = xl["R"].values
d1 = xl["Dperp"].values
r2 = xl["R.3"].values
d2 = xl["Dperp.3"].values
f1 = interp1d(r1, d1)
f2 = interp1d(r2, d2)

fig, ax = plt.subplots(figsize=(5, 4))
ax.axvline(2.3125, color="k", linestyle="--", lw=2, zorder=5)
ax.fill_between(r1, d1, f2(r1), color="tab:red", alpha=0.5, zorder=15)
ax.plot(r1, d1, color="tab:red", lw=2, zorder=15)
ax.plot(r1, f2(r1), color="tab:red", lw=2, zorder=15)
ax.set_yscale("log")
ax.grid(alpha=0.3, which="both")
ax.set_ylabel(r"$\mathdefault{D^{eff}_{r}=\Gamma_r\ /\ (dn/dr)\ \ (m^2/s)}$", fontsize=12)
ax.set_xlabel("Radial (m)", fontsize=12)
fig.tight_layout()
fig.show()