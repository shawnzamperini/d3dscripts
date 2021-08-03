import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from numpy.polynomial import Polynomial

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Input parameters.
charge      = 15
poly_order  = 7
#ring        = 17

# 167277-inj-006 inputs.
#ring = 17
stag_region = [15, 50]; play = 5; ring = 17; root_num = 0 # For 006 ring 17.
#stag_region = [15, 50]; play = 10
#stag_region = [20, 30]; play = 7; ring = 17; root_num = 0 # For 006d ring 17.

# 167247-inj-034 inputs (Dperp scan).
#stag_region = [30, 50]; play = 10; ring = 30; root_num = 2

# What case.
path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc"
#path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034a.nc"
op = oedge_plots.OedgePlots(path)

# Also load all the cases bc why not.
ops = []
cases = ["006", "006a", "006b", "006c", "006d", "006e"]
for case in cases:
    path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-{:}.nc".format(case)
    ops.append(oedge_plots.OedgePlots(path))

# Load the along ring goodies.
s, vz = op.along_ring(ring, "VELAVG", charge=charge, plot_it=False)
s, nz = op.along_ring(ring, "DDLIMS", charge=charge, plot_it=False)
gz = nz * vz

# Polynomial fits around the stagnation point. First fit is to fine tune the
# stagnation point.
stag1 = np.logical_and(s >= stag_region[0], s <= stag_region[1])
p = Polynomial(poly_order).fit(s[stag1], vz[stag1], poly_order)
roots = p.roots()
root = roots[np.isreal(roots)].real[root_num]
print("Stagnation at {:.3f}".format(root))

# Second fits are the real ones.
stag2 = np.logical_and(s >= root-play, s <= root+play)
pv = Polynomial(poly_order).fit(s[stag2], vz[stag2], poly_order)
pn = Polynomial(poly_order).fit(s[stag2], nz[stag2], poly_order)
#pg = Polynomial(poly_order).fit(s[stag2], gz[stag2], poly_order)
pg = pn * pv

# Calculate the criteria for accumulation.
gz_d3 = pg.deriv(3)(root)
vz_d3 = pv.deriv(3)(root)
vz_d1 = pv.deriv(3)(root)
criteria = (gz_d3 - pn(root) * vz_d3) / (3 * vz_d1)
print("Criteria = {:.3e}".format(criteria))

# Derivative functions for the alternative criteria.
dgz = pg.deriv(1)
dvz = pv.deriv(1)
dgz2 = pg.deriv(2)
dvz2 = pv.deriv(2)

# Maximum values.
vz_max = np.abs(vz[stag2]).max()
nz_max = np.abs(nz[stag2]).max()
gz_max = np.abs(gz[stag2]).max()

# Some plots to vizualize.
x = np.linspace(s[stag2].min(), s[stag2].max(), 100)
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.axhline(0.0, color="k", linestyle="-")
ax1.axvline(root, color="k", linestyle="-")
ax1.plot(s[stag2], vz[stag2]/vz_max, color="r", label="vz")
ax1.plot(s[stag2], nz[stag2]/nz_max, color="g", label="nz")
ax1.plot(s[stag2], gz[stag2]/gz_max, color="b", label="gz")
ax1.plot(x, pv(x)/vz_max, color="r", linestyle="--")
ax1.plot(x, pn(x)/nz_max, color="g", linestyle="--")
ax1.plot(x, pg(x)/gz_max, color="b", linestyle="--")
crit2 = (dgz(x) - pg(x) / pv(x) * dvz(x)) / pv(x)
crit3 = 1 / pv(x) * dgz2(x) - pg(x) / (pv**2)(x) * dvz2(x) + 2 * pg(x) / (pv**3)(x) * dvz(x) * dvz(x)
ax1.plot(x, crit2/crit2.max(), color="cyan", linestyle="-")
ax1.plot(x, crit3/crit3.max(), color="magenta", linestyle="-")
ax1.set_xlabel("Distance from inner target (m)")
ax1.set_ylabel("Normalized values")
ax1.legend()

stag_points = [42.6, 42.0, 38.7, 29.9, 26.1, 22.6]
for i in range(0, len(ops)):
    s, vz = ops[i].along_ring(ring, "VELAVG", charge=charge, plot_it=False)
    s, nz = ops[i].along_ring(ring, "DDLIMS", charge=charge+1, plot_it=False)
    #stag_point = s[vz[np.where(vz!=0)[0]]]
    ax2.axvline(stag_points[i], color="C{}".format(i), linestyle="--")
    #ax2.axhline(0, color="k")
    ax2.plot(s, nz, color="C{}".format(i))

fig.tight_layout()
fig.show()
