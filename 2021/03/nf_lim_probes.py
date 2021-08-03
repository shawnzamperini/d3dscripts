# This script well generate plots of the deposition profiles and their comparisons
# to RBS and LAMS data.
import sys
sys.path.append('/mnt/c/Users/Shawn/Documents/GitHub/d3dscripts/2020/05')
sys.path.append('/mnt/c/Users/Shawn/Documents/GitHub/d3dscripts/2020/08')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

fav_shift = -1.0

# Load in A8 data and pull out relevant arrays.
import compare_a8_3dlim
a8_lim_itfotf_x = compare_a8_3dlim.lim_com_x
a8_lim_itfotf_y = compare_a8_3dlim.lim_itfotf
a8_lam_itfotf_x = compare_a8_3dlim.lams_com_x
a8_lam_itfotf_y = compare_a8_3dlim.lams_itfotf
a8_lam_itf_x = compare_a8_3dlim.lams_itf_x
a8_lam_itf_y = compare_a8_3dlim.lams_itf_y
a8_lam_otf_x = compare_a8_3dlim.lams_otf_x
a8_lam_otf_y = compare_a8_3dlim.lams_otf_y
a8_rbs_itf_x = compare_a8_3dlim.rbs_itf_x
a8_rbs_itf_y = compare_a8_3dlim.rbs_itf_y
a8_rbs_otf_x = compare_a8_3dlim.rbs_otf_x
a8_rbs_otf_y = compare_a8_3dlim.rbs_otf_y
a8_lim_itf_x = compare_a8_3dlim.lim_itf_x
a8_lim_itf_y = compare_a8_3dlim.lim_itf_y
a8_lim_otf_x = compare_a8_3dlim.lim_otf_x
a8_lim_otf_y = compare_a8_3dlim.lim_otf_y

# Load in A15 data and pull out relevant arrays.
import compare_a15_3dlim
a15_lim_itfotf_x = compare_a15_3dlim.lim_com_x
a15_lim_itfotf_y = compare_a15_3dlim.lim_itfotf
a15_lam_itfotf_x = compare_a15_3dlim.lams_com_x
a15_lam_itfotf_y = compare_a15_3dlim.lams_itfotf
a15_lam_itf_x = compare_a15_3dlim.lams_itf_x
a15_lam_itf_y = compare_a15_3dlim.lams_itf_y
a15_lam_otf_x = compare_a15_3dlim.lams_otf_x
a15_lam_otf_y = compare_a15_3dlim.lams_otf_y
a15_rbs_itf_x = compare_a15_3dlim.rbs_itf_x
a15_rbs_itf_y = compare_a15_3dlim.rbs_itf_y
a15_rbs_otf_x = compare_a15_3dlim.rbs_otf_x
a15_rbs_otf_y = compare_a15_3dlim.rbs_otf_y
a15_lim_itf_x = compare_a15_3dlim.lim_itf_x_all
a15_lim_itf_y = compare_a15_3dlim.lim_itf_y_all
a15_lim_otf_x = compare_a15_3dlim.lim_otf_x_all
a15_lim_otf_y = compare_a15_3dlim.lim_otf_y_all

# Renormalize all the data for each probe.
a8_lim_itf_y = a8_lim_itf_y / np.max((a8_lim_itf_y.max(), a8_lim_otf_y.max()))
a8_lim_otf_y = a8_lim_otf_y / np.max((a8_lim_itf_y.max(), a8_lim_otf_y.max()))
a8_lam_itf_y = a8_lam_itf_y / np.max((a8_lam_itf_y.max(), a8_lam_otf_y.max()))
a8_lam_otf_y = a8_lam_otf_y / np.max((a8_lam_itf_y.max(), a8_lam_otf_y.max()))
a8_rbs_itf_y = a8_rbs_itf_y / np.max((a8_rbs_itf_y.max(), a8_rbs_otf_y.max()))
a8_rbs_otf_y = a8_rbs_otf_y / np.max((a8_rbs_itf_y.max(), a8_rbs_otf_y.max()))
a15_lim_itf_y = a15_lim_itf_y / np.max((a15_lim_itf_y.max(), a15_lim_otf_y.max()))
a15_lim_otf_y = a15_lim_otf_y / np.max((a15_lim_itf_y.max(), a15_lim_otf_y.max()))
a15_lam_itf_y = a15_lam_itf_y / np.max((a15_lam_itf_y.max(), a15_lam_otf_y.max()))
a15_lam_otf_y = a15_lam_otf_y / np.max((a15_lam_itf_y.max(), a15_lam_otf_y.max()))
a15_rbs_itf_y = a15_rbs_itf_y / np.max((a15_rbs_itf_y.max(), a15_rbs_otf_y.max()))
a15_rbs_otf_y = a15_rbs_otf_y / np.max((a15_rbs_itf_y.max(), a15_rbs_otf_y.max()))

# Plot a grid of both sides of the probes.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 16
lw = 3
band_alpha = 0.6
ms = 12
xlims = [5.5, 15.5]

# Century Gothic is missing a unicode minus sign so we gotta create the tick
# labels ourselves.
ytick_labels = [r"10$^{-\mathdefault{2}}$", r"10$^{-\mathdefault{1}}$",
  r"10$^{-\mathdefault{0}}$"]

fig = plt.figure(figsize=(9, 6))
ax  = fig.add_subplot(111)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# Background figure so that we can share labels.
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

ax1.axvline(11.5, color="k", linestyle="--", lw=lw, zorder=1)
ax1.axvline(7.5, color="k", linestyle="--", lw=lw, zorder=2)
ax1.fill_between(a8_lam_itf_x, (a8_lam_itf_y+a8_lam_itf_y*0.2),
  (a8_lam_itf_y-a8_lam_itf_y*0.2), alpha=band_alpha, color=colors[2], zorder=2,
  label="LAMS")
ax1.plot(a8_lim_itf_x, a8_lim_itf_y, color='k', lw=lw+2)
ax1.plot(a8_lim_itf_x, a8_lim_itf_y, color=colors[2], lw=lw, label="3DLIM")
ax1.plot(a8_rbs_itf_x, a8_rbs_itf_y, '*', color=colors[2], ms=ms, mec='k', label="RBS")
ax1.set_yscale("log")
ax1.set_xlim(xlims)
ax1.set_ylim([0.01, 2])
ax1.grid()
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.tick_params(axis="x", which="both", labelbottom=False)
ax1.set_title("Unfavorable", fontsize=fontsize)
ax1.tick_params(axis="y", which="both", labelsize=12)
#ax1.set_xticks(np.arange(5, 16, 2))
ax1.set_yticks([0.01, 0.1, 1])
ax1.set_yticklabels(ytick_labels)
#ax1.legend()
ax1.text(0.05, 0.9, "ITF", transform=ax1.transAxes, fontsize=fontsize)
ax1.text(0.04, 0.05, "UDL", transform=ax1.transAxes, fontsize=fontsize)
ax1.text(0.3, 0.05, "UBL", transform=ax1.transAxes, fontsize=fontsize)
ax1.text(0.8, 0.13, "OWL", transform=ax1.transAxes, fontsize=fontsize)

ax3.fill_between(a8_lam_otf_x, (a8_lam_otf_y+a8_lam_otf_y*0.2),
  (a8_lam_otf_y-a8_lam_otf_y*0.2), alpha=band_alpha, color=colors[2])
ax3.plot(a8_lim_otf_x, a8_lim_otf_y, color='k', lw=lw+2)
ax3.plot(a8_lim_otf_x, a8_lim_otf_y, color=colors[2], lw=lw)
ax3.plot(a8_rbs_otf_x, a8_rbs_otf_y, '*', color=colors[2], ms=ms, mec='k')
ax3.set_yscale("log")
ax3.set_xlim(xlims)
ax3.set_ylim([0.01, 2])
ax3.grid()
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)
ax3.set_xlabel(r"R-$\mathdefault{R_{sep}}$ OMP (cm)", fontsize=fontsize)
ax.set_ylabel("W Deposition (Normalized)\n", fontsize=fontsize)
ax3.tick_params(axis="both", which="both", labelsize=12)
#ax3.set_xticks(np.arange(5, 16, 2))
ax3.set_yticks([0.01, 0.1, 1])
ax3.set_yticklabels(ytick_labels)
ax3.text(0.05, 0.9, "OTF", transform=ax3.transAxes, fontsize=fontsize)

ax2.axvline(11.5, color="k", linestyle="--", lw=lw, zorder=1)
ax2.axvline(7.5, color="k", linestyle="--", lw=lw, zorder=2)
#ax2.axvline(8.5, color="k", linestyle="--", lw=lw, zorder=1)
#ax2.axvline(11.5, color="k", linestyle="--", lw=lw, zorder=1)
ax2.fill_between(a15_lam_itf_x+fav_shift, (a15_lam_itf_y+a15_lam_itf_y*0.2),
  (a15_lam_itf_y-a15_lam_itf_y*0.2), alpha=band_alpha, color=colors[3],
  label="LAMS")
ax2.plot(a15_lim_itf_x+fav_shift, a15_lim_itf_y, color='k', lw=lw+2)
ax2.plot(a15_lim_itf_x+fav_shift, a15_lim_itf_y, color=colors[3], lw=lw, label="3DLIM")
ax2.plot(a15_rbs_itf_x+fav_shift, a15_rbs_itf_y, '*', color=colors[3], ms=ms, mec='k', label="RBS")
ax2.set_yscale("log")
ax2.set_xlim(xlims)
ax2.set_ylim([0.01, 2])
ax2.grid()
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.tick_params(axis="x", which="both", labelbottom=False)
ax2.set_title("Favorable", fontsize=fontsize)
ax2.tick_params(axis="x", which="both", labelbottom=False)
ax2.tick_params(axis="y", which="both", labelleft=False)
#ax2.set_xticks(np.arange(5, 16, 2))
ax2.legend()
ax2.text(0.05, 0.9, "ITF", transform=ax2.transAxes, fontsize=fontsize)
#ax2.text(0.05, 0.9, "ITF", transform=ax2.transAxes, fontsize=fontsize)
ax2.text(0.04, 0.05, "UDL", transform=ax2.transAxes, fontsize=fontsize)
ax2.text(0.3, 0.05, "UBL", transform=ax2.transAxes, fontsize=fontsize)
ax2.text(0.63, 0.05, "OWL", transform=ax2.transAxes, fontsize=fontsize)

ax4.fill_between(a15_lam_otf_x+fav_shift, (a15_lam_otf_y+a15_lam_otf_y*0.2),
  (a15_lam_otf_y-a15_lam_otf_y*0.2), alpha=band_alpha, color=colors[3])
ax4.plot(a15_lim_otf_x+fav_shift, a15_lim_otf_y, color='k', lw=lw+2)
ax4.plot(a15_lim_otf_x+fav_shift, a15_lim_otf_y, color=colors[3], lw=lw)
ax4.plot(a15_rbs_otf_x+fav_shift, a15_rbs_otf_y, '*', color=colors[3], ms=ms, mec='k')
ax4.set_yscale("log")
ax4.set_xlim(xlims)
ax4.set_ylim([0.01, 2])
ax4.grid()
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)
ax4.set_xlabel(r"R-$\mathdefault{R_{sep}}$ OMP (cm)", fontsize=fontsize)
ax4.tick_params(axis="x", which="both", labelsize=13)
ax4.tick_params(axis="y", which="both", labelleft=False)
#ax4.set_xticks(np.arange(5, 16, 2))
ax4.text(0.05, 0.9, "OTF", transform=ax4.transAxes, fontsize=fontsize)

fig.tight_layout()
fig.show()

# Now a plot of the ITF/OTF ratio for each probe.
fig, ax = plt.subplots(figsize=(7, 4))

ax.axvline(11.5, color="k", linestyle="--", lw=1, zorder=1)
ax.axvline(7.5, color="k", linestyle="--", lw=1, zorder=2)
ax.grid(zorder=5)
ax.axhline(1.0, color="k", linestyle="--", lw=1)
ax.fill_between(a8_lam_itfotf_x, a8_lam_itfotf_y+a8_lam_itfotf_y*0.1,
  a8_lam_itfotf_y-a8_lam_itfotf_y*0.1, alpha=band_alpha, color=colors[2], lw=lw,
  zorder=6)
ax.fill_between(a15_lam_itfotf_x+fav_shift, a15_lam_itfotf_y+a15_lam_itfotf_y*0.1,
  a15_lam_itfotf_y-a15_lam_itfotf_y*0.1, alpha=band_alpha, color=colors[3], lw=lw, zorder=7)
ax.plot(a8_lim_itfotf_x[10:],  a8_lim_itfotf_y[10:], '-', color="k", lw=lw+2, zorder=8)
ax.plot(a8_lim_itfotf_x[10:],  a8_lim_itfotf_y[10:], '-', color=colors[2], lw=lw, zorder=9)
ax.plot(a15_lim_itfotf_x[10:]+fav_shift,  a15_lim_itfotf_y[10:], '-', color="k", lw=lw+2, zorder=10)
ax.plot(a15_lim_itfotf_x[10:]+fav_shift,  a15_lim_itfotf_y[10:], '-', color=colors[3], lw=lw, zorder=11)
ax.set_ylim([0, 4])
ax.set_xlim(xlims)
ax.set_xlabel(r"R-$\mathdefault{R_{sep}}$ OMP (cm)", fontsize=fontsize)
ax.set_ylabel("ITF/OTF", fontsize=fontsize)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_yticks(np.arange(0, 5))
ax.text(0.05, 0.9, "UDL", transform=ax.transAxes, fontsize=fontsize)
ax.text(0.38, 0.9, "UBL", transform=ax.transAxes, fontsize=fontsize)
ax.text(0.8, 0.6, "OWL", transform=ax.transAxes, fontsize=fontsize)

thick_black = Line2D([0], [0], color='k', lw=lw+2)
a8_line = Line2D([0], [0], color=colors[2], lw=lw)
a15_line = Line2D([0], [0], color=colors[3], lw=lw)
custom_lines = ((thick_black, a8_line), (thick_black, a15_line))
ax.legend(custom_lines, ["Unfavorable","Favorable"], fontsize=fontsize)

ax.annotate("3DLIM", (8.33, 1.85), xytext=(9.1, 1.5),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize,
  bbox=dict(facecolor="white", ec="none"))
ax.annotate("LAMS", (8.31, 2.66), xytext=(7.5, 3.2),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize,
  bbox=dict(facecolor="white", ec="none"))
ax.annotate("LAMS", (8.15, 0.62), xytext=(7.0, 0.2),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize,
  bbox=dict(facecolor="white", ec="none"))
ax.annotate("3DLIM", (9.57, 0.33), xytext=(9.5, 0.7),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize,
  bbox=dict(facecolor="white", ec="none"))

fig.tight_layout()
fig.show()
