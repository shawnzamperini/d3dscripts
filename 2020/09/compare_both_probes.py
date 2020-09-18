import sys
sys.path.append('/mnt/c/Users/Shawn/Documents/GitHub/d3dscripts/2020/05')
sys.path.append('/mnt/c/Users/Shawn/Documents/GitHub/d3dscripts/2020/08')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


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
fig, ax = plt.subplots(figsize=(10,6))
ax.axhline(1, linestyle='--', color='k', lw=2)
ax.plot(a8_lam_itfotf_x,  a8_lam_itfotf_y, color="C4", lw=lw, label=r"$\mathrm{Bx\nabla B\ Away}$")
#ax.plot(a8_lim_itfotf_x,  a8_lim_itfotf_y, '--', color="C4", lw=lw)

# Normal dashed line.
#ax.plot(a8_lim_itfotf_x[10:],  a8_lim_itfotf_y[10:], '--', color="C4", lw=lw)

# Solid line with outline.
ax.plot(a8_lim_itfotf_x[10:],  a8_lim_itfotf_y[10:], '-', color="k", lw=lw+2)
ax.plot(a8_lim_itfotf_x[10:],  a8_lim_itfotf_y[10:], '-', color="C4", lw=lw)

ax.plot(a15_lam_itfotf_x, a15_lam_itfotf_y, color="C5", lw=lw, label=r"$\mathrm{Bx\nabla B\ Towards}$")

# Normal dashed line.
ax.plot(a15_lim_itfotf_x, a15_lim_itfotf_y, '--', color="C5", lw=lw)

# Solid line with outline.
ax.plot(a15_lim_itfotf_x, a15_lim_itfotf_y, '-', color="k", lw=lw+2)
ax.plot(a15_lim_itfotf_x, a15_lim_itfotf_y, '-', color="C5", lw=lw)

ax.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
ax.set_ylabel("ITF/OTF", fontsize=fontsize)
ax.tick_params(which='both', labelsize=labelsize)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticks(np.arange(0, 10, 2))
ax.set_yticks(np.arange(0, 6, 1))
ax.legend(custom_lines, ['3DLIM', 'LAMS'], fontsize=fontsize)
ax.set_xlim([0, 7])
ax.set_ylim([0, 4])
fig.tight_layout()
fig.show()
