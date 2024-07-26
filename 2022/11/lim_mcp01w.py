import LimPlots
import matplotlib.pyplot as plt
import pandas as pd

# absfac = 7.124e12 * 3.5  # ABSFAC * duration of flattop
absfac = 2e16
ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/190423/190423-mcp01w-002.nc"
lp = LimPlots.LimPlots(ncpath)

cent = lp.centerline(showplot=False)
lim_itfx = cent["x2"]
lim_itfy = cent["y2"] * absfac
lim_otfx = cent["x1"]
lim_otfy = cent["y1"] * absfac

# The LAMS data.
itf_path = "/Users/zamperini/My Drive/Research/Data/cp_data/MCPL01W.csv"
otf_path = "/Users/zamperini/My Drive/Research/Data/cp_data/MCPR01W.csv"
itf = pd.read_csv(itf_path)
otf = pd.read_csv(otf_path)
lam_itfx = itf["Distance along probe (cm)"].values / 100
lam_itfy = itf["W areal density (1e15cm-2)"].values * 1e15 * 10000
lam_otfx = otf["Distance along probe (cm)"].values / 100
lam_otfy = otf["W areal density (1e15cm-2)"].values * 1e15 * 10000

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

ax1.plot(lam_itfx, lam_itfy, color="tab:red", label="LAMS", lw=2)
ax1.plot(lim_itfx, lim_itfy, color="k", lw=3)
ax1.plot(lim_itfx, lim_itfy, color="tab:red", label="3DLIM", lw=2)
ax1.set_ylabel("W Areal Density (m-3)")
ax1.set_xlabel("Distance along probe (m)")
ax1.set_title("ITF")

ax2.plot(lam_otfx, lam_otfy, color="tab:purple", label="LAMS", lw=2)
ax2.plot(lim_otfx, lim_otfy, color="k", lw=3)
ax2.plot(lim_otfx, lim_otfy, color="tab:purple", label="3DLIM", lw=2)
ax2.set_xlabel("Distance along probe (m)")
ax2.set_title("OTF")

fig.tight_layout()
fig.show()
plt.show()
