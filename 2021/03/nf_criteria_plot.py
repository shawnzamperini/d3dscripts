# This script will plot the density and the criteria for accumulation.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Input parameters.
charge      = 15
#poly_order  = 11

path1 = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc"
path2 = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006d.nc"
op1 = oedge_plots.OedgePlots(path1)
op2 = oedge_plots.OedgePlots(path2)
#stag_region = [15, 50]; play = 5; ring = 17; root_num = 0

# Set diffusion coefficients.
dperp1 = 0.3
dperp2 = 0.3

ops = [op1, op2]
xs = []; ss = []; nzs_raw = []; nzs_fit = []; vzs_raw = []; vzs_fit = []
crits = []; gzs_fit = []; dgzs_fit = []; crits2 = []; dnz2s_fit = []; crits3 = []
crit3_play = 0; zeros = [42.5, 25.8]
for i in range(0, len(ops)):

    # Case specific parameters.
    if ops[i] == op1:
        stag_region = [15, 50]; play = 5; ring = 17; root_num = 0; poly_order = 7
        #stag_region = [37, 50]; play = 5; ring = 17; root_num = 0; poly_order = 3
    elif ops[i] == op2:
        stag_region = [20, 30]; play = 8; ring = 17; root_num = 0; poly_order = 7

    # Load the along ring goodies.
    s, vz = ops[i].along_ring(ring, "VELAVG", charge=charge, plot_it=False)
    s, nz = ops[i].along_ring(ring, "DDLIMS", charge=charge, plot_it=False)
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
    #weights = np.abs(root / s[stag2])
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

    # For the plots.
    x = np.linspace(s[stag2].min(), s[stag2].max(), 100)

    # The second erivative of nz as derived.
    dnz2 = 1/pv(x)*dgz2(x) - 2/(pv**2)(x)*dgz(x)*dvz(x) - pg(x)/(pv**2)(x)*dvz2(x) + 2*pg(x)/(pv**3)(x)*(dvz**2)(x)

    # Our criteria via the second derivative. Replace avlues near the singularity with a nan.
    crit3 = dvz2(x) / pv(x) - dgz2(x) / pg(x)
    #crit3 = dgz2(x) / dvz2(x)
    crit3[np.logical_and(x<zeros[i]+crit3_play, x>zeros[i]-crit3_play)] = np.nan

    crit4 = dgz2(x) / dvz2(x)
    #fig, ax = plt.subplots()
    #ax.plot(x, pn(x), color="k")
    #ax.plot(x, crit4, color="r")
    #fig.tight_layout()
    #fig.show()

    ss.append(s[stag2])
    xs.append(x)
    nzs_raw.append(nz[stag2]/nz_max)
    nzs_fit.append(pn(x)/nz_max)
    vzs_raw.append(vz[stag2]/vz_max)
    vzs_fit.append(pv(x)/vz_max)
    gzs_fit.append(pg(x)/gz_max)
    dgzs_fit.append(dgz(x)/np.abs(dgz(x)).max())
    crit = (dgz(x) - pg(x) / pv(x) * dvz(x)) / pv(x)  # Just dn/ds is all this is.
    crits.append(crit/crit.max())
    crit2 = 2 * pg(x) / (pv**3)(x) * (dvz**2)(x) - 2 / (pv*2)(x) * dgz(x) * dvz(x)
    crits2.append(crit2/np.abs(crit2).max())
    dnz2s_fit.append(dnz2/np.abs(dnz2.max()))
    crits3.append(crit3/np.nanmax(np.abs(crit3)))

# Distance to separatrix at OMP.
mid_dist = ops[0].nc.variables["MIDIST"][1][ring]
mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(mid_dist*100)

# Plotting.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 5
combine = False

# The plots.
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 8))

ax1.axhline(0, color="k")
ax1.plot(xs[0], crits[0],   color=colors[2], label=r"$\mathdefault{dn_Z}$/ds", lw=lw)
ax1.plot(ss[0], vzs_raw[0], color=colors[3], alpha=0.4, lw=lw)
ax1.plot(xs[0], vzs_fit[0], color=colors[3], label=r"$\mathdefault{v_Z}$", lw=lw)
ax1.plot(ss[0], nzs_raw[0], color=colors[1], alpha=0.4, lw=lw)
ax1.plot(xs[0], nzs_fit[0], color=colors[1], label=r"$\mathdefault{n_Z}$", lw=lw)
ax1.plot(xs[0], crits3[0],  color=colors[4], label="Crit3", lw=lw)
ax1.legend(fontsize=12)
ax1.grid(zorder=1)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.set_ylabel("Normalized values", fontsize=fontsize)
#ax1.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax1.text(0.05, 0.1, mid_str, transform=ax1.transAxes, fontsize=fontsize, bbox=dict(color="white"))
ax1.text(0.05, 0.95, "a) W15+ no additional flow", transform=ax1.transAxes, fontsize=fontsize, bbox=dict(color="white"))
ax1.set_ylim([-1, 1.2])

ax2.axhline(0, color="k")
ax2.plot(xs[1], crits[1],   color=colors[2], label=r"$\mathdefault{dn_Z}$/ds", lw=lw)
ax2.plot(ss[1], vzs_raw[1], color=colors[3], alpha=0.4, lw=lw)
ax2.plot(xs[1], vzs_fit[1], color=colors[3], label=r"$\mathdefault{v_Z}$", lw=lw)
ax2.plot(ss[1], nzs_raw[1], color=colors[1], alpha=0.4, lw=lw)
ax2.plot(xs[1], nzs_fit[1], color=colors[1], label=r"$\mathdefault{n_Z}$", lw=lw)
ax2.plot(xs[1], crits3[1],  color=colors[4], label="Crit3", lw=lw)
ax2.grid(zorder=1)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.set_ylabel("Normalized values", fontsize=fontsize)
ax2.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax2.text(0.05, 0.95, "b) W15+ M = 0.4 additional flow", transform=ax2.transAxes, fontsize=fontsize, bbox=dict(color="white"))
#ax2.set_ylim([-0.6, 1.2])
ax2.set_ylim([-1.0, 1.2])

fig.tight_layout()
fig.show()
