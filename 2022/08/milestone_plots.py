import oedge_plots
import matplotlib.pyplot as plt
import LimPlots
import pandas as pd
from scipy.signal import savgol_filter, medfilt


ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-sput-001.nc"
op = oedge_plots.OedgePlots(ncpath)

s, nz = op.along_ring(13, "DDLIMS", charge="all", plot_it=False)

rmrsomp = float(op.nc["MIDIST"][1][13-1])

fig, ax = plt.subplots(figsize=(5,4))
ax.plot(s, nz, lw=3, color="tab:red")
ax.set_title("R-Rsep @ OMP = {:.1f} mm".format(rmrsomp*1000), fontsize=14)
ax.set_xlabel("Distance from W rings (m)", fontsize=14)
ax.set_ylabel("W Density (arb)", fontsize=14)
#ax.set_yscale("log")
#ax.set_ylim(1e13, 1e15)
ax.set_ylim([0, 0.5e15])
ax.grid(which="both")

fig.tight_layout()
fig.show()

# Aaron's simulations.
ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/190423/190423-test-8.i.nc"
#ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/190423/190423-test-8.i-2.pdv0.1.nc" # Mach = 0.1
#ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/190423/190423-test-8.i-2.pdv0.2.nc" # Mach = 0.2
lp = LimPlots.LimPlots(ncpath)
cent = lp.centerline(showplot=False)
rorigin = ((2.295 + 0.002) - 2.2339) * 100
#rshift = 0.6
rshift = 0
limitfx = rorigin + cent["x2"] * 100 + rshift
limitfy = cent["y2"]
limotfx = rorigin + cent["x1"] * 100 + rshift
limotfy = cent["y1"]

itfpath = "/Users/zamperini/Documents/d3d_work/cpdata/MCPL01W.csv"
otfpath = "/Users/zamperini/Documents/d3d_work/cpdata/MCPR01W.csv"
itf = pd.read_csv(itfpath)
otf = pd.read_csv(otfpath)
itfx = itf["Distance along probe (cm)"] + rorigin
itfy = itf["W areal density (1e15cm-2)"]
otfx = otf["Distance along probe (cm)"] + rorigin
otfy = otf["W areal density (1e15cm-2)"]

itfy = medfilt(itfy, 21)
otfy = medfilt(otfy, 21)
limitfy = savgol_filter(limitfy, 15, 2)
limotfy = savgol_filter(limotfy, 15, 2)

c1 = "tab:orange"
c2 = "tab:blue"
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4), sharex=True, sharey=True)

idx = 5
ax1.plot(itfx[idx:], itfy[idx:], color=c2, lw=3, label="ITF")
ax2.plot(otfx[idx:], otfy[idx:], color=c1, lw=3, label="OTF")

ax11 = ax1.twinx()
ax22 = ax2.twinx()
ax11.plot(limitfx[:-5], limitfy[:-5], color="k", lw=4)
ax11.plot(limitfx[:-5], limitfy[:-5], color=c2, lw=3)
ax22.plot(limotfx[5:], limotfy[5:], color="k", lw=4)
ax22.plot(limotfx[5:], limotfy[5:], color=c1, lw=3)

ax1.legend()
ax2.legend()
fig.supxlabel("R-Rsep (cm)", fontsize=14)
ax1.set_ylabel("W Areal Density (1e15 cm-2)", fontsize=14)
ax22.set_ylabel("Simulated Deposition (arb)", fontsize=14)
ax1.set_xlim([6, 12])
#ax1.set_ylim([0, 0.025])
#ax11.set_ylim([0, 50])
#ax22.set_ylim([0, 50])
ax1.set_ylim([0.001, 0.03])
ax11.set_ylim([1, 50])
ax22.set_ylim([1, 50])
ax1.set_yscale("log")
ax2.set_yscale("log")
ax11.set_yscale("log")
ax22.set_yscale("log")
ax1.grid(alpha=0.3, which="both")
ax2.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.show()
