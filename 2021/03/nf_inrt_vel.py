# Point of this script is to demonstrate how including an inertial velocity
# obtains significantly better agreement with the W velocity.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from numpy.polynomial import Polynomial

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Some constants.
#ring = 40
mz_amu = 183.84
md_amu = 2.01
mz = mz_amu * 1.66e-27  # Mass of W in kg
md = md_amu * 1.66e-27   # Mass of D in kg
col_log = 15
#charges = np.arange(1, 31)
charges = [15]
ring = 17
#ring = 30

# Just do a basic case.
fav_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006d.nc"
#fav_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034a.nc"
fav = oedge_plots.OedgePlots(fav_path)

# Average vz weighted by their densities.
def weighted_vz(op, ring):
    all_vz = []; all_nz = []; all_tau = []
    for charge in charges:
        s, vz = op.along_ring(ring, "VELAVG", charge=charge, plot_it=False)
        s, nz = op.along_ring(ring, "DDLIMS", charge=charge, plot_it=False)
        s, ti = op.along_ring(ring, "KTIBS", plot_it=False)
        s, ne = op.along_ring(ring, "KNBS",  plot_it=False)
        tau_s = 1.47E13 * mz_amu * ti * np.sqrt(ti / md_amu) / \
                ((1 + md_amu / mz_amu) * ne * np.power(charge, 2) * col_log)
        all_vz.append(vz)
        all_nz.append(nz)
        all_tau.append(tau_s)
    all_vz = np.array(all_vz)
    all_nz = np.array(all_nz)
    all_tau = np.array(all_tau)

    # Remove all charges where no data exists.
    #mask1 = all_nz.sum(axis=1)!=0.0
    #all_vz = all_vz[mask1]
    #all_nz = all_nz[mask1]
    #all_tau = all_tau[mask1]

    # Remove all locations where no data exists.
    mask2 = all_nz.sum(axis=0)!=0.0
    all_vz = all_vz[:, mask2]
    all_nz = all_nz[:, mask2]
    all_tau = all_tau[:, mask2]
    s = s[mask2]

    # Weighted averages.
    #avg_vz = np.average(all_vz, weights=all_nz, axis=0)
    #avg_nz = np.average(all_nz, weights=all_nz, axis=0)
    #avg_tau = np.average(all_tau, weights=all_nz, axis=0)
    avg_vz = all_vz[0]
    avg_nz = all_nz[0]
    avg_tau = all_tau[0]

    # Term used in the inertial velocity.
    gamma = savgol_filter(avg_vz, 21, 2) * savgol_filter(avg_nz, 21, 2)
    tmp = savgol_filter(np.power(savgol_filter(avg_vz, 21, 2), 2) * savgol_filter(avg_nz, 21, 2), 21, 2)
    vint = avg_tau / avg_nz * np.gradient(savgol_filter(gamma * avg_vz, 21, 2), s)
    #vint = avg_tau / avg_nz * np.gradient(np.power(avg_vz, 2) * avg_nz, s)

    return {"s":s, "avg_vz":avg_vz, "avg_nz":avg_nz, "mask":mask2, "vint":vint}

# Need to use the mask on the other variables so it can line up with the data from
# the weighting process.
_, vi = fav.along_ring(ring, "Velocity", plot_it=False)
_, ti = fav.along_ring(ring, "KTIBS", plot_it=False)
_, ne = fav.along_ring(ring, "KNBS", plot_it=False)

# Now grab the inertial velocity and some others.
weighted = weighted_vz(fav, ring)
s_w    = weighted["s"]
avg_vz = weighted["avg_vz"]
avg_nz = weighted["avg_nz"]
avg_gamma = avg_vz * avg_nz
mask   = weighted["mask"]
vint   = weighted["vint"]
vi_w = vi[mask]
ti_w = ti[mask]
ne_w = ne[mask]

# Just set charge = 1 since they will cancel out in the vti calculation.
tau_s = 1.47E13 * mz_amu * ti_w * np.sqrt(ti_w / md_amu) / \
        ((1 + md_amu / mz_amu) * ne_w * np.power(1, 2) * col_log)
_, fig = fav.along_ring(ring, "FIG", charge=1, plot_it=False)
fig_w = fig[mask]
vti = fig_w * tau_s / mz

# Distance to separatrix at OMP.
mid_dist = fav.nc.variables["MIDIST"][1][ring]
mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(mid_dist*100)

# Plot it up.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 5
comb_str = r"$\mathdefault{v_i+v_{Ti}}$" + "-" + r"$\mathdefault{v_{inertia}}$"

vz_max = np.max((np.abs(avg_vz.min()), np.abs(avg_vz.max())))
nz_max = np.max((np.abs(avg_nz.min()), np.abs(avg_nz.max())))
gamma_max = np.max((np.abs(avg_gamma.min()), np.abs(avg_gamma.max())))

fig, ax1 = plt.subplots(figsize=(5.5, 4))

ax1.axhline(0.0, linestyle="-", color="k")
#ax1.plot(s_w, vi_w + vti, label=r"$\mathdefault{v_i+v_{Ti}}$", color=colors[3], lw=lw)
#ax1.plot(s_w, -vint, linestyle="--", color=colors[3], lw=lw)
#ax1.plot(s_w, vi_w + vti - vint, label=comb_str, color=colors[2], lw=lw)
ax1.plot(s_w, avg_gamma/gamma_max, label="gamma", color=colors[2], lw=lw)
ax1.plot(s_w, avg_vz/vz_max, label=r"$\mathdefault{v_Z}$", color=colors[0], lw=lw)
ax1.plot(s_w, avg_nz/nz_max, label=r"$\mathdefault{n_Z}$", color=colors[1], lw=lw)
ax1.legend(fontsize=14, loc="upper right")
ax1.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
#ax1.set_ylabel("W8+ Velocity (m/s)", fontsize=fontsize)
ax1.set_ylabel("Normalized Values", fontsize=fontsize)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
#ax1.set_ylim([-25000, 25000])
ax1.grid()
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.text(0.05, 0.95, mid_str, fontsize=fontsize, bbox=dict(color="white"), transform=ax1.transAxes)

fig.tight_layout()
fig.show()

# Fit to polynomials. Limit to just the region around the stagnation point.
stag_start = 20; stag_end = 30
stag_region = np.logical_and(s_w > stag_start, s_w < stag_end)
poly_order = 5
p = Polynomial(poly_order).fit(s_w[stag_region], avg_gamma[stag_region], poly_order)

# Find the root (stagnation point) and refit.
roots = p.roots()
root = roots[np.isreal(roots)].real[0]
print("Stagnation at {}".format(root))
stag_region = np.logical_and(s_w > root-3, s_w < root+3)
pg = Polynomial(poly_order).fit(s_w[stag_region], avg_gamma[stag_region], poly_order)
pv = Polynomial(poly_order).fit(s_w[stag_region], avg_vz[stag_region], poly_order)
pn = Polynomial(poly_order).fit(s_w[stag_region], avg_nz[stag_region], poly_order)

# Derivatives for my criteria.
gamma_d3 = pg.deriv(3)(root)
vz_d3 = pv.deriv(3)(root)
vz_d1 = pv.deriv(3)(root)
#print("Gamma``` = {}".format(gamma_d3))
criteria = (gamma_d3 - pn(root) * vz_d3) / 3 * vz_d1
print("Criteria = {:.3f}".format(criteria))

fig, ax1 = plt.subplots()

ax1.axhline(0.0, color="k", linestyle="--")
ax1.plot(s_w[stag_region], avg_gamma[stag_region])
ax1.plot(s_w[stag_region], pg(s_w[stag_region]))
ax1.set_xlabel("Distance from inner target (m)")
ax1.set_ylabel("Gamma")

fig.tight_layout()
fig.show()
