import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import splrep, BSpline
import numpy as np
from scipy.signal import medfilt
from matplotlib import lines
import matplotlib.patches as mpatches


fname = "/Users/zamperini/github/d3dscripts/2024/01/imp_inner_data.pickle"
with open(fname, "rb") as f:
    data = pickle.load(f)

def spline_fit(x, y, s, log=False, npoints=100):
    if log:
        tck = splrep(x, np.log(y), s=s)
        xnew = np.linspace(x.min(), x.max(), npoints)
        return xnew, np.exp(BSpline(*tck)(xnew))
    else:
        tck = splrep(x, y, s=s)
        xnew = np.linspace(x.min(), x.max(), npoints)
        return xnew, BSpline(*tck)(xnew)

vrhist_exb_r = data["vrhist_exb_r"]
vrhist_pol_r = data["vrhist_pol_r"]
rhist = data["rhist"]
time = data["time"]
t_spline, pol_spline = spline_fit(time, medfilt(vrhist_pol_r, 11), 1e9, npoints=500)
t_spline, exb_spline = spline_fit(time, vrhist_exb_r, 1e8, npoints=500)
t_spline, rhist_spline = spline_fit(time, rhist, 1e-3, npoints=500)


fig, (ax1, ax11) = plt.subplots(1, 2, figsize=(9, 4), sharex=True)
#ax11 = ax1.twinx()
exb_scale = 1
ax1.axhline(0, color="k")
# ax1.plot(time*1e6, vrhist_exb_r * exb_scale, color="tab:red", alpha=0.3)
ax1.plot(t_spline*1e6, exb_spline * exb_scale, color="k", lw=4)
ax1.plot(t_spline*1e6, exb_spline * exb_scale, color="tab:red", lw=3)
# ax1.plot(time*1e6, vrhist_pol_r, color="tab:purple", alpha=0.3)
ax1.plot(t_spline*1e6, pol_spline, color="k", lw=4)
ax1.plot(t_spline*1e6, pol_spline, color="tab:purple", lw=3)
ax1.set_ylim([-2e4, 1e4])
# ax11.plot(time*1e6, rhist, color="k", alpha=0.3)
ax11.plot(t_spline*1e6, rhist_spline, color="k", lw=4)
# ax11.axhline((rmin + bound_buffer) - gk_rsep, color="k", linestyle="--")
ax1.set_xlabel("Time (us)", fontsize=14)
ax11.set_xlabel("Time (us)", fontsize=14)
ax1.set_ylabel("Radial Velocity (m/s)", fontsize=14)
ax11.set_ylabel(r"R-$\mathdefault{R_{sep}}$ (m)", fontsize=14)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax11.tick_params(axis='both', which='major', labelsize=12)

def add_arrow(x_tail, y_tail, x_head, y_head, color, mut):
    dx = x_head - x_tail
    dy = y_head - y_tail
    arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (x_head, y_head),
                                 mutation_scale=mut, color=color)
    ax1.add_patch(arrow)

add_arrow(2.5, 5000, 6, 5000, "k", 20)
add_arrow(3, 1500, 0, 1500, "tab:red", 20)
add_arrow(3, -4000, 0, -4000, "tab:purple", 20)

# The legend.
exb_line1 = lines.Line2D([], [], linewidth=3, linestyle="-", color="tab:red")
exb_line2 = lines.Line2D([], [], linewidth=4, linestyle="-", color="k")
pol_line1 = lines.Line2D([], [], linewidth=3, linestyle="-", color="tab:purple")
pol_line2 = lines.Line2D([], [], linewidth=4, linestyle="-", color="k")
ax1.legend([(exb_line2, exb_line1), (pol_line2, pol_line1)],
           [r"$\mathdefault{v_r^{ExB}}$", r"$\mathdefault{v_r^{pol}}$"], fontsize=14, loc="lower left")

fig.tight_layout()
fig.show()
