# Make a plot of the SiC yields for my SiC 4 pager.
import matplotlib.pyplot as plt
import sys
import numpy as np

sys.path.append("../../2022/08")
import EngelhardtModel
em = EngelhardtModel.EngelhardtModel()
em.load_mm_model()

ep = np.geomspace(25, 1000, 1000)
ec = np.geomspace(5, 1000, 1000)

lw = 3
fontsize = 12
fig, ax1 = plt.subplots(figsize=(4, 3.5))

ax1.plot(ep, em.Y_D_C(ep), color="tab:red", label="D->C", lw=lw)
ax1.plot(ec, em.Y_D_C_ch(ec), color="tab:red", linestyle="--", label="D->C (ch)", lw=lw)
ax1.plot(ep, em.Y_D_SiC_C(ep), color="tab:green", label="D->SiC", lw=lw)
ax1.plot(ec, em.Y_D_SiC_Cch(ec), color="tab:green", linestyle="--", label="D->SiC (ch)", lw=lw)
ax1.plot(ep, em.Y_D_SiC_Si(ep), color="tab:purple", label="D->SiC,Si", lw=lw)

ax1.annotate("Chemical", (15, em.Y_D_C_ch(15)), (5, 0.07), rotation=0, fontsize=13, backgroundcolor="w",
             # bbox={"facecolor":"w", "edgecolor":"w", "boxstyle":"round"},
             arrowprops={"arrowstyle":"-", "color":"k", "linestyle":"--", "lw":2})
ax1.annotate("Physical", (400, em.Y_D_C(400)), (150, 0.3), rotation=0, fontsize=13, backgroundcolor="w",
             # bbox={"facecolor":"w", "edgecolor":"w", "boxstyle":"round"},
             arrowprops={"arrowstyle":"-", "color":"k", "linestyle":"-", "lw":2})

# ax1.text(670, em.Y_D_C(670), "C", fontsize=11, bbox={"facecolor":"tab:red", "edgecolor":"k"})
# ax1.text(670, em.Y_D_C(670), "C", fontsize=11, bbox={"facecolor":"tab:red", "edgecolor":"k"})
# ax1.text(670, em.Y_D_C(670), "C", fontsize=11, bbox={"facecolor":"tab:red", "edgecolor":"k"})

ax1.set_ylabel("Yield", fontsize=fontsize)
ax1.set_xlabel("Impact Energy (eV)", fontsize=fontsize)
ax1.set_title("Deuterium Impact", fontsize=fontsize)
ax1.set_ylim([1e-5, 1.0])
ax1.grid(which="both", alpha=0.3)
#ax1.legend()
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.text(1.02, 0.04, "Si (SiC)", transform=ax1.transAxes, rotation=90, color="tab:purple", fontsize=fontsize)
ax1.text(1.02, 0.32, "C (SiC)", transform=ax1.transAxes, rotation=90, color="tab:green", fontsize=fontsize)
ax1.text(1.02, 0.60, "C (Graphite)", transform=ax1.transAxes, rotation=90, color="tab:red", fontsize=fontsize)
fig.tight_layout()
fig.show()