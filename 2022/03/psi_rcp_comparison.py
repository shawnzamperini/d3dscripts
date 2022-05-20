import pandas as pd
import matplotlib.pyplot as plt


mimes = True

rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"
if mimes:
    rcp267 = pd.read_excel(rcp_path, sheet_name="MP184267_1")
    rcp527 = pd.read_excel(rcp_path, sheet_name="MP184527_1")
    xcol = "R (cm)"
else:
    rcp267 = pd.read_excel(rcp_path, sheet_name="XP184267_1")
    rcp527 = pd.read_excel(rcp_path, sheet_name="XP184527_1")
    xcol = "Z (cm)"

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, figsize=(8, 8))

x267 = rcp267[xcol] / 100
x527 = rcp527[xcol] / 100

c267 = "tab:red"
c527 = "tab:purple"

cut267 = 2.354
cut527 = 2.344

ylim = [-0.01, 0.105]
ax1.axhline(0, color="k", linestyle="--")
ax2.axhline(0, color="k", linestyle="--")
ax1.plot(x267, rcp267["Isat (A)"], color=c267)
ax2.plot(x527, rcp527["Isat (A)"], color=c527)
ax1.set_ylabel("Isat (A)")
ax1.set_ylim(ylim)
ax2.set_ylim(ylim)
ax1.set_title("#184267", color=c267)
ax2.set_title("#184527", color=c527)
if mimes:
    ax1.axvline(cut267, color="k", linestyle=":")
    ax2.axvline(cut527, color="k", linestyle=":")

# 184267: +M = ITF
# 184527: +M = OTF, so we multiply by -1
ylim = [-0.4, 0.75]
ax3.axhline(0, color="k", linestyle="--")
ax4.axhline(0, color="k", linestyle="--")
ax3.plot(x267, rcp267["Mach"], color=c267)
ax4.plot(x527, -rcp527["Mach"], color=c527)
ax3.set_ylabel("Mach")
ax3.set_ylim(ylim)
ax4.set_ylim(ylim)
if mimes:
    ax3.axvline(cut267, color="k", linestyle=":")
    ax4.axvline(cut527, color="k", linestyle=":")
    ax3.text(2.3, 0.6, "Towards Inner")
    ax3.text(2.3, -0.2, "Towards Outer")
    ax4.text(2.3, 0.6, "Towards Inner")
    ax4.text(2.3, -0.2, "Towards Outer")

ylim = [-1, 25]
ax5.axhline(0, color="k", linestyle="--")
ax6.axhline(0, color="k", linestyle="--")
ax5.plot(x267, rcp267["Te (eV)"], color=c267)
ax6.plot(x527, rcp527["Te (eV)"], color=c527)
ax5.set_ylabel("Te (eV)")
ax5.set_ylim(ylim)
ax6.set_ylim(ylim)
if mimes:
    ax5.axvline(cut267, color="k", linestyle=":")
    ax6.axvline(cut527, color="k", linestyle=":")

ylim = [-2, 8]
ax7.axhline(0, color="k", linestyle="--")
ax8.axhline(0, color="k", linestyle="--")
ax7.plot(x267, rcp267["ne (1e18 m-3)"], color=c267)
ax8.plot(x527, rcp527["ne (1e18 m-3)"], color=c527)
ax7.set_ylabel("ne (1e18 m-3)")
ax7.set_ylim(ylim)
ax8.set_ylim(ylim)
if mimes:
    ax7.axvline(cut267, color="k", linestyle=":")
    ax8.axvline(cut527, color="k", linestyle=":")


fig.tight_layout()
fig.show()
