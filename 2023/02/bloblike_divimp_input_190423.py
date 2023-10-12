import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np


ca1 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_190411_1.tab", delimiter="\t")
ca2 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_190411_2.tab", delimiter="\t")

# Some additional calculations.
dt = 0.005
fblob1 = ca1["Npeaks"] / dt
fblob2 = ca2["Npeaks"] / dt
vrbar1 = ca1["Vr(m/s)"] * ca1["T_blob(e-6s)"] * 1e-6 * fblob1
vrbar2 = ca2["Vr(m/s)"] * ca2["T_blob(e-6s)"] * 1e-6 * fblob2

# Create a distribution assuming normally distributed values.
mean_vr = pd.concat((ca1, ca2))["Vr(m/s)"].mean()
std = mean_vr / 2
vrs = np.linspace(0, mean_vr * 3, 150)
rv = norm()
pdf = norm.pdf(vrs, loc=mean_vr, scale=std)

# Print out for copy/paste into DIVIMP option
for i in range(0, len(vrs)):
    print("{:.2f} {:.2e}".format(vrs[i], pdf[i]))

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(vrs, pdf)
ax.set_xlabel("vr (m/s)")
ax.set_ylabel("prob")
fig.tight_layout()
fig.show()

# Plot summary.
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 4), sharex=True)

ax1.axvline(0.0, color="k", linestyle="--")
ax1.plot(ca1["R-Rsep(cm)"], ca1["Vr(m/s)"], label="#1", color="tab:red")
ax1.plot(ca2["R-Rsep(cm)"], ca2["Vr(m/s)"], label="#2", color="tab:purple")
x = [0, 10]
y = np.full(len(x), mean_vr)
yerr = np.full(len(x), std)
ax1.fill_between(x, y-yerr, y+yerr, alpha=0.3, color="k")
ax1.axhline(mean_vr, color="k")
ax1.set_ylabel("vr (m/s)")
ax1.set_xlim([-0.5, 9])
ax1.legend()

ax2.axvline(0.0, color="k", linestyle="--")
ax2.plot(ca1["R-Rsep(cm)"], fblob1, label="#1", color="tab:red")
ax2.plot(ca2["R-Rsep(cm)"], fblob2, label="#2", color="tab:purple")
ax2.set_ylabel("fblob (Hz)")

ax3.axvline(0.0, color="k", linestyle="--")
ax3.plot(ca1["R-Rsep(cm)"], ca1["D_rad(cm)"], label="#1", color="tab:red")
ax3.plot(ca2["R-Rsep(cm)"], ca2["D_rad(cm)"], label="#2", color="tab:purple")
ax3.set_ylabel("Drad (cm)")
ax3.set_xlabel("R-Rsep (cm)")

ax4.axvline(0.0, color="k", linestyle="--")
ax4.plot(ca1["R-Rsep(cm)"], vrbar1, label="#1", color="tab:red")
ax4.plot(ca2["R-Rsep(cm)"], vrbar2, label="#2", color="tab:purple")
ax4.set_ylabel("vrbar (m/s)")
ax4.set_xlabel("R-Rsep (cm)")

fig.tight_layout()
fig.show()
