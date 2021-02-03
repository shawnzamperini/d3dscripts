import matplotlib.pyplot as plt
import lim_plots as lim
from scipy.signal import savgol_filter


ncpatha = "/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-005a.nc"
ncpathb = "/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-005b.nc"
ncpathc = "/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-005c.nc"

lpa = lim.LimPlots(ncpatha)
lpb = lim.LimPlots(ncpathb)
lpc = lim.LimPlots(ncpathc)
da = lpa.centerline(show_plot=False)
db = lpb.centerline(show_plot=False)
dc = lpc.centerline(show_plot=False)

x = da["otf_x"] * 100
ya = da["otf_y"]
yb = db["otf_y"]
yc = dc["otf_y"]

# Smooth the data some.
ya = savgol_filter(ya, 25, 2)
yb = savgol_filter(yb, 25, 2)
yc = savgol_filter(yc, 25, 2)

fig, ax = plt.subplots(figsize=(8,5))
ax.plot(x, ya, label="7 m", lw=4, color="coral")
ax.plot(x, yb, label="5 m", lw=4, color="springgreen")
ax.plot(x, yc, label="3 m", lw=4, color="orchid")
ax.set_xlabel("Distance along probe (cm)", fontsize=16)
ax.set_ylabel("Deposition (arbitrary)", fontsize=16)
ax.tick_params(which='both', labelsize=12)
ax.set_xlim([3, 13])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_yscale("log")
ax.legend(fontsize=16)
fig.tight_layout()
fig.show()
