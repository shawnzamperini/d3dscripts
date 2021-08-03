import sys
sys.path.append('/mnt/c/Users/Shawn/Documents/GitHub/d3dscripts/2020/05')
sys.path.append('/mnt/c/Users/Shawn/Documents/GitHub/d3dscripts/2020/08')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple


band = True
band_alpha = 0.6
log = False
to_omp = False

# Load in A8 data.
import compare_a8_3dlim
a8_lim_itfotf_x = compare_a8_3dlim.lim_com_x
a8_lim_itfotf_y = compare_a8_3dlim.lim_itfotf
a8_lam_itfotf_x = compare_a8_3dlim.lams_com_x
a8_lam_itfotf_y = compare_a8_3dlim.lams_itfotf

# Load in A15 data.
import compare_a15_3dlim
a15_lim_itfotf_x = compare_a15_3dlim.lim_com_x
a15_lim_itfotf_y = compare_a15_3dlim.lim_itfotf
a15_lam_itfotf_x = compare_a15_3dlim.lams_com_x
a15_lam_itfotf_y = compare_a15_3dlim.lams_itfotf

# Needed for the legend.
custom_lines = [Line2D([0], [0], color='k', lw=3, linestyle='--'), Line2D([0], [0], color='k', lw=3)]

# Create a plot comparing the two.
plt.rcParams['font.family'] = 'sans-serif'
plt.style.use("tableau-colorblind10")
lw = 5
fontsize  = 16
labelsize = 11
a8_color = "orange"
a15_color = "lightskyblue"

thick_black = Line2D([0], [0], color='k', lw=lw+2)
a8_line = Line2D([0], [0], color=a8_color, lw=lw)
a15_line = Line2D([0], [0], color=a15_color, lw=lw)
custom_lines = ((thick_black, a8_line), (thick_black, a15_line), a8_line, a15_line)

if to_omp:
    xlabel = "R-Rsep OMP (cm)"
    xbounds = [6.5, 14.5]
    xlabels = np.arange(7, 17, 2)
else:
    xlabel = "Distance along probe (cm)"
    xbounds = [0, 7]
    xlabels = np.arange(0, 10, 2)

fig, ax = plt.subplots(figsize=(10,6))
ax.axhline(1, linestyle='--', color='k', lw=2)
if band:
    ax.fill_between(a8_lam_itfotf_x, a8_lam_itfotf_y+a8_lam_itfotf_y*0.1, a8_lam_itfotf_y-a8_lam_itfotf_y*0.1, alpha=band_alpha, color=a8_color)
else:
    ax.plot(a8_lam_itfotf_x,  a8_lam_itfotf_y, color=a8_color, lw=lw, label=r"$\mathrm{Bx\nabla B\ Away}$")
#ax.plot(a8_lim_itfotf_x,  a8_lim_itfotf_y, '--', color="C4", lw=lw)

# Normal dashed line.
#ax.plot(a8_lim_itfotf_x[10:],  a8_lim_itfotf_y[10:], '--', color="C4", lw=lw)

# Solid line with outline.
ax.plot(a8_lim_itfotf_x[10:],  a8_lim_itfotf_y[10:], '-', color="k", lw=lw+2)
ax.plot(a8_lim_itfotf_x[10:],  a8_lim_itfotf_y[10:], '-', color=a8_color, lw=lw)

if band:
    ax.fill_between(a15_lam_itfotf_x, a15_lam_itfotf_y+a15_lam_itfotf_y*0.1, a15_lam_itfotf_y-a15_lam_itfotf_y*0.1, alpha=band_alpha, color=a15_color)
else:
    ax.plot(a15_lam_itfotf_x, a15_lam_itfotf_y, color=a15_color, lw=lw, label=r"$\mathrm{Bx\nabla B\ Towards}$")

# Normal dashed line.
ax.plot(a15_lim_itfotf_x, a15_lim_itfotf_y, '--', color=a15_color, lw=lw)

# Solid line with outline.
ax.plot(a15_lim_itfotf_x, a15_lim_itfotf_y, '-', color="k", lw=lw+2)
ax.plot(a15_lim_itfotf_x, a15_lim_itfotf_y, '-', color=a15_color, lw=lw)

#ax.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
ax.set_xlabel(xlabel, fontsize=fontsize)
ax.set_ylabel("ITF/OTF", fontsize=fontsize)
ax.tick_params(which='both', labelsize=labelsize)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_xticks(np.arange(0, 10, 2))
ax.set_xticks(xlabels)
ax.set_yticks(np.arange(0, 6, 1))
ax.legend(custom_lines, ["","","", ""], fontsize=fontsize)
#ax.set_xlim([0, 7])
ax.set_xlim(xbounds)
ax.set_ylim([0, 4])
if log:
    ax.set_yscale("log")
    ax.set_ylim([0.1, 5])
fig.tight_layout()
fig.show()
