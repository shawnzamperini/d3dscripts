import pandas as pd
import matplotlib.pyplot as plt



# Don't forget to mount g drive first
#  sudo mount -t drvfs G: /mnt/g
path = "/mnt/g/My Drive/Research/Documents/2026/flan_friction_force_test.xlsx"
df = pd.read_excel(path, skiprows=7, sheet_name="therm_samp")

# Load data into arrays. First set is calculating s with mean_g, second set is 
# with inst_g.
mean_t = df["t"].values
mean_vx = df["vX"].values
mean_vx_std = df["vX std"].values
mean_anly = df["analytic"].values
inst_t = df["t.1"].values
inst_vx = df["vX.1"].values
inst_vx_std = df["vX std.1"].values
inst_anly = df["analytic.1"].values

fontsize = 16
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.fill_between(mean_t*1e6, mean_vx-mean_vx_std, mean_vx+mean_vx_std, 
	color="tab:red", alpha=0.3)
ax1.plot(mean_t*1e6, mean_vx, color="tab:red", lw=4, label="Mean")

ax1.fill_between(inst_t*1e6, inst_vx-inst_vx_std, inst_vx+inst_vx_std, 
	color="tab:purple", alpha=0.3)
ax1.plot(inst_t*1e6, inst_vx, color="tab:purple", lw=4, label="Instant")

ax1.plot(mean_t*1e6, mean_anly, color="k", lw=3, linestyle="--")

ax1.set_xlabel(r"Time ($\mathdefault{\mu s}$)", fontsize=fontsize)
ax1.set_ylabel(r"v (m/s)", fontsize=fontsize)
ax1.tick_params(axis='both', labelsize=fontsize-2)
ax1.legend(fontsize=fontsize)

fig.tight_layout()
fig.show()

