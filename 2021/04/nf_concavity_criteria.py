# This script goes from about the X-point to the crown area and estimates the
# radial gamma at each location by fitting a polynomial to the radial values.
# It uses the weighted /(by the charge densities) averages of vz in these calculations.
# It then fits polynomials to the parallel velocty and gamma, and plug them
# all into the criteria for negative concavity.
import oedge_plots
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial
from scipy.signal import savgol_filter
from matplotlib import ticker
from scipy.optimize import curve_fit

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Order for polynomials.
poly_order = 7
window = 15; sf_poly = 3

# SOL rings to calculate the gradient across.
ring1 = 19; ring2 = 70
num_rings = ring2 - ring1

# Knots that go from about the X-point to the crown region.
knot1 = 65; knot2 = 120  # Fine tune these.
num_knots = knot2 - knot1

charges = np.arange(0, 30)  # First 30 charge states.

def exp_fit(x, a, b, c):
    return a *  np.exp(b * x) + c

dperps = [0.3, 0.4, 0.6, 1.0, 5.0]
dp03_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034a2.nc"  # Running more particles on dipsy.
dp06_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034g.nc"
dp1_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034b.nc"
dp5_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034c.nc"
dp10_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034d.nc"

dp03 = oedge_plots.OedgePlots(dp03_path)
dp06 = oedge_plots.OedgePlots(dp06_path)
dp1 = oedge_plots.OedgePlots(dp1_path)
dp5 = oedge_plots.OedgePlots(dp5_path)
dp10 = oedge_plots.OedgePlots(dp10_path)
ops = [dp03, dp06, dp1, dp5, dp10]

# Same sized array that for each column of the same knot are the distances
# from the first/closest knot.
dists = np.zeros((num_rings, num_knots))
rs = dp03.nc.variables["RS"][ring1-1:ring2-1, knot1-1:knot2-1].data
zs = dp03.nc.variables["ZS"][ring1-1:ring2-1, knot1-1:knot2-1].data
for i in range(0, num_rings):
    for j in range(0, num_knots):
        p1 = (rs[0, j], zs[0, j])  # The point nearest to the separatrix.
        p2 = (rs[i, j], zs[i, j])
        dists[i, j] = np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

# Lists to hold everything.
all_vzs = []; all_nzs = []; all_crits = []; all_nzs_poly_par = []
all_nzs_poly_rad = []; all_vzs_poly_par = []; all_vzs_lin = []
all_vzs_slopes = []; all_slope_crits = []; all_nzs_vals = []
all_nzsd2_sav_rad = []

#for o in range(1, 2):
for o in range(0, len(ops)):
    op = ops[o]
    print(op)

    # First fill an array of weighted averages of vz. Dimensions
    # are (charge, ring, knot).
    print("  Loading along ring data...")
    ss  = np.zeros((num_rings, num_knots))
    vzs = np.zeros((len(charges), num_rings, num_knots))
    nzs = np.zeros((len(charges), num_rings, num_knots))
    for i in range(0, num_rings):
        ring = ring1 + i
        print("   Ring: {}".format(ring))
        for j in range(0, len(charges)):
            charge = charges[j]
            s, vz = op.along_ring(ring, "VELAVG", charge=charge, plot_it=False, remove_zeros=False)
            s, nz = op.along_ring(ring, "DDLIMS", charge=charge, plot_it=False, remove_zeros=False)
            vzs[j, i] = vz[knot1:knot2]
            nzs[j, i] = nz[knot1:knot2]
        ss[i] = s[knot1:knot2]

    print("  Calculating weighted averages and parallel polynomials...")
    vzs_weighted = np.zeros((num_rings, num_knots))
    nzs_summed   = np.zeros((num_rings, num_knots))
    nzs_poly_par = []; vzs_poly_par = []; vzs_lin = []
    for i in range(0, num_rings):
        for j in range(0, num_knots):
            nzs_summed[i, j] = nzs[:, i, j].sum()
            if nzs_summed[i, j] != 0.0:
                vzs_weighted[i, j] = np.average(vzs[:, i, j], weights=nzs[:, i, j])

        # Fit polynomials to the parallel density and velocity. Radial done later.
        nz_poly = Polynomial(poly_order).fit(ss[i], nzs_summed[i],   poly_order)
        vz_poly = Polynomial(poly_order).fit(ss[i], vzs_weighted[i], poly_order)
        nzs_poly_par.append(nz_poly)
        vzs_poly_par.append(vz_poly)

        # Linear fit to the velocity.
        vz_poly = Polynomial(1).fit(ss[i], vzs_weighted[i], 1)
        vzs_lin.append(vz_poly)

    nzs_poly_par = np.array(nzs_poly_par)
    vzs_poly_par = np.array(vzs_poly_par)

    print("  Calculating radial polynomials...")
    nzs_poly_rad = []; vzs_poly_rad = [];
    nzs_d2_exp_rad = np.zeros((num_rings, num_knots));
    nzs_d2_sav_rad = np.zeros((num_rings, num_knots))
    for j in range(0, num_knots):
        rad_dist = dists[:, j]
        rad_nzs  = nzs_summed[:, j]
        #rad_vzs  = vzs_weighted[:, j]
        nz_poly = Polynomial(poly_order).fit(rad_dist, rad_nzs, poly_order)
        popt, pcov = curve_fit(exp_fit, rad_dist, rad_nzs, maxfev=5000)
        grad1 = np.gradient(savgol_filter(rad_nzs, 11, 2), rad_dist)
        grad2 = np.gradient(savgol_filter(grad1, 11, 2), rad_dist)
        #vz_poly = Polynomial(poly_order).fit(rad_dist, rad_vzs, poly_order)
        nzs_poly_rad.append(nz_poly)
        nzs_d2_sav_rad[:, j] = grad2

        # The second derivative of the exponential fit.
        d2 = popt[0] * popt[1]**2 * np.exp(popt[1] * rad_dist)
        nzs_d2_exp_rad[:, j] = d2

        #vzs_poly_rad.append(vz_poly)
    nzs_poly_rad = np.array(nzs_poly_rad)
    #vzs_poly_rad = np.array(vzs_poly_rad)

    # Using the polynomials calculate the respective gammas.
    print("  Calculating gammas...")
    dg_poly_par = []; dg_poly_rad = []
    dgs_par = np.zeros((num_rings, num_knots))
    dgs_rad = np.zeros((num_rings, num_knots))
    dgs_rad_exp = np.zeros((num_rings, num_knots))
    dgs_rad_sav = np.zeros((num_rings, num_knots))
    dg_crit = np.zeros((num_rings, num_knots))
    dg_crit_exp = np.zeros((num_rings, num_knots))
    dg_crit_sav = np.zeros((num_rings, num_knots))
    for i in range(0, num_rings):
        pn_par = nzs_poly_par[i]
        pv_par = vzs_poly_par[i]
        dpn_par = pn_par.deriv(1)
        dpv_par = pv_par.deriv(1)
        dpg_par = pv_par * dpn_par + pn_par * dpv_par
        for j in range(0, num_knots):
            pn_rad = nzs_poly_rad[j]
            dpn2_rad = pn_rad.deriv(2)
            dpg_rad = - dperps[o] * dpn2_rad
            #dg_crit[i, j] = dpg_par(ss[i, j]) + dpg_rad(dists[i, j])
            dgs_par[i, j] = dpg_par(ss[i, j])
            dgs_rad[i, j] = dpg_rad(dists[i, j])
            dgs_rad_exp[i, j] = -dperps[o] * nzs_d2_exp_rad[i, j]
            dgs_rad_sav[i, j] = -dperps[o] * nzs_d2_sav_rad[i, j]
            dg_crit_exp[i, j] = dgs_par[i, j] + dgs_rad_exp[i, j]
            dg_crit_sav[i, j] = dgs_par[i, j] + dgs_rad_sav[i, j]

        # Calculate the criteria here so we can savgol smooth dpg_rad.
        #dgs_rad[i] = savgol_filter(dgs_rad[i], window, sf_poly)
        dg_crit[i] = dgs_par[i] + savgol_filter(dgs_rad[i], window, sf_poly)

    # Plot to see how the parallel polynomials look.
    #every = 25
    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))
    #ax1.plot(ss[::every], nzs_summed[::every])
    #ax2.plot(ss[::every], vzs_weighted[::every])
    #for i in np.arange(0, len(nzs_poly_par), every):
    #    ax1.plot(ss[::every], nzs_poly_par[i](ss[::every]))
    #    ax2.plot(ss[::every], vzs_poly_par[i](ss[::every]))
    #ax1.set_xlabel("Distance from inner target (m)")
    #ax1.set_ylabel("Density (arbitrary)")
    #ax2.set_xlabel("Distance from inner target (m)")
    #ax2.set_ylabel("Velocity (m/s)")
    #fig.tight_layout()
    #fig.show()

    # Append to lists for saving.
    all_vzs.append(vzs_weighted)
    all_nzs.append(nzs_summed)
    all_crits.append(dg_crit)
    all_nzs_poly_par.append(nzs_poly_par)
    all_nzs_poly_rad.append(nzs_poly_rad)
    all_vzs_poly_par.append(vzs_poly_par)
    all_vzs_lin.append(vzs_lin)
    #all_nzs_exp_rad.append(nzs_exp_rad)
    all_nzsd2_sav_rad.append(nzs_d2_sav_rad)

    # Plot to see how the criteria look compared to the density.
    plot_ring = 30 - ring1 - 1; plot_knot = 90 - knot1 - 1
    fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))
    #ax.plot(ss[plot_ring], dg_crit[plot_ring]/np.abs(dg_crit[plot_ring]).max(), color="tab:red", label="criteria")
    #ax.plot(ss[plot_ring], dg_crit_exp[plot_ring]/np.abs(dg_crit_exp[plot_ring]).max(), color="tab:red", label="criteria")
    ax.plot(ss[plot_ring], savgol_filter(dg_crit_sav[plot_ring], window, sf_poly)/np.abs(savgol_filter(dg_crit_sav[plot_ring], window, sf_poly)).max(), color="tab:red", label="criteria")
    ax.plot(ss[plot_ring], nzs_summed[plot_ring]/np.abs(nzs_summed[plot_ring]).max(), color="tab:purple", alpha=0.5)
    poly_y = nzs_poly_par[plot_ring](ss[plot_ring])
    ax.plot(ss[plot_ring], poly_y/np.abs(nzs_summed[plot_ring]).max(), color="tab:purple", label="nz")
    ax.plot(ss[plot_ring], vzs_weighted[plot_ring]/np.abs(vzs_weighted[plot_ring]).max(), color="tab:green", alpha=0.5)
    poly_y = vzs_poly_par[plot_ring](ss[plot_ring])
    ax.plot(ss[plot_ring], poly_y/np.abs(vzs_weighted[plot_ring]).max(), color="tab:green", label="vz")
    ax.set_xlabel("Distance from inner target (m)")
    #ax.set_ylabel("Criteria")
    ax.legend()

    vzs_slopes = []; nzs_vals = []; slope_crits = []
    for i in range(0, num_rings):
        l = vzs_lin[i]

        # Calculate the zero point.
        b = l.coef[0]
        m = l.coef[1]
        root = -b / m
        vzs_slopes.append(m)
        stag_idx = np.abs(vzs_lin[i](ss[i])).argmin()

        # Is the density above or below average.
        #nzs_val = nzs_poly_par[i](root) / nzs_poly_par[i](ss[i]).mean()
        nzs_val = nzs_summed[i][stag_idx] / nzs_summed[i].mean()
        nzs_vals.append(nzs_val)

        #slope_crit = dgs_rad[i][stag_idx] / nzs_poly_par[i](root)
        #slope_crit = -dgs_rad[i][stag_idx] / nzs_summed[i][stag_idx]
        slope_crit = -dgs_rad_exp[i][stag_idx] / nzs_summed[i][stag_idx]
        slope_crits.append(slope_crit)
    vzs_slopes = np.array(vzs_slopes)
    slope_crits = np.array(slope_crits)
    nzs_vals = np.array(nzs_vals)
    all_vzs_slopes.append(vzs_slopes)
    all_slope_crits.append(slope_crits)
    all_nzs_vals.append(nzs_vals)

    line_x = np.linspace(np.min(vzs_slopes), np.max(vzs_slopes), 3)
    #ax2.plot(line_x, line_x, color="k", linestyle="--")
    #ax2.scatter(vzs_slopes, slope_crits)
    ax2.scatter(vzs_slopes-slope_crits, nzs_vals)
    ax2.set_xlabel("vz_slope-slope_crit")
    ax2.set_ylabel("nz/nz_avg")

    """
    rad = dists[ring1-1:ring2-1, plot_knot]
    rad_nz = nzs_summed[ring1-1:ring2-1, plot_knot]
    ax2.plot(rad, rad_nz/np.abs(rad_nz).max(), color="tab:purple", alpha=0.5)
    ax2.plot(rad, nzs_poly_rad[plot_knot](rad)/np.abs(rad_nz).max(), color="tab:purple", label="nz")
    poly_y = nzs_poly_rad[plot_knot].deriv(1)(rad)/np.abs(nzs_poly_rad[plot_knot].deriv(1)(rad)).max()
    ax2.plot(rad, poly_y, color="tab:purple", label="dnz/ds")
    poly_y = nzs_poly_rad[plot_knot].deriv(2)(rad)/np.abs(nzs_poly_rad[plot_knot].deriv(2)(rad)).max()
    ax2.plot(rad, poly_y, color="tab:purple", label="d2nz/ds2")
    ax2.set_xlabel("Distance from separatrix (m)")
    ax2.legend()
    """

    ax3.plot(ss[plot_ring], dgs_par[plot_ring], color="b", label="dg/ds")
    ax3.plot(ss[plot_ring], dgs_rad[plot_ring], color="r", alpha=0.5)
    ax3.plot(ss[plot_ring], savgol_filter(dgs_rad[plot_ring], window, sf_poly), color="r", label="dg/dr")
    ax3.plot(ss[plot_ring], dgs_rad_sav[plot_ring], color="g", alpha=0.5)
    ax3.plot(ss[plot_ring], savgol_filter(dgs_rad_sav[plot_ring], window, sf_poly), color="g", label="dg/dr_sav")
    #ax3.plot(ss[plot_ring], nzs_poly_par[plot_ring].deriv(2)(ss[plot_ring]), color="g", label="d2nz/ds2")
    ax3.set_xlabel("Distance from inner target (m)")
    ax3.set_ylabel("Gamma derivative")
    ax3.legend()

    fig.tight_layout()
    fig.show()

    # Heatmaps of the density and criteria.
#    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))
#    nzs_norm = np.zeros(nzs_summed.shape)
#    for i in range(0, num_rings):
#        nzs_norm[i] = nzs_summed[i]/nzs_summed[i].min()
#    gt = dg_crit >= 0; lt = dg_crit < 0
#    dg_norm = np.zeros(dg_crit.shape)
#    dg_norm = dg_crit / np.abs(dg_crit).max()
#    #cs1 = ax1.contourf(nzs_norm, locator=ticker.LogLocator())
#    cs1 = ax1.contourf(nzs_norm)
#    levels = np.linspace(-1, 1, 10)
#    cs2 = ax2.contourf(dg_norm, cmap="coolwarm", levels=levels)
#    cbar1 = fig.colorbar(cs1)
#    cbar2 = fig.colorbar(cs2)
#    fig.tight_layout()
#    fig.show()

# Plots to see how the weighted averages look.
#fig, ax = plt.subplots(figsize=(7, 7))
#ax.plot(ss[::2], vzs_weighted[::2])
#ax.contourf(dists)
#ax.set_xlabel("Distance from inner target (m)")
#ax.set_ylabel("Velocity (m/s)")
#fig.tight_layout()
#fig.show()

# Plot to see how the radial velocity looks.
#fig, ax = plt.subplots(figsize=(7, 7))
#ax.plot(dists[::2], nzs_summed[::2])
#ax.contourf(nzs_summed)
#ax.set_xlabel("Distance from separatrix")
#ax.set_ylabel("Density (arbitrary)")
#fig.tight_layout()
#fig.show()

# Could do a plot of the crits vs. the densities / min density.
fig, ax1 = plt.subplots(figsize=(8, 5))

for i in range(0, len(ops)):
    vzs_slopes = all_vzs_slopes[i]
    slopes_crits = all_slope_crits[i]
    nzs_vals = all_nzs_vals[i]
    ax1.scatter(vzs_slopes - slopes_crits, nzs_vals, label=i)

ax1.legend()
ax1.set_xlabel("c - crit")
ax1.set_ylabel("nz/nz_avg")
fig.tight_layout()
fig.show()

def lin_fit(x, m, b):
    return m * x + b

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 5

# Plot of the parallel derivative of gammaz for each case.
fig, ax1 = plt.subplots(figsize=(5.5, 4))

xlim1 = 36; xlim2 = 42
ax1.axvline(39, color="k", linestyle="--", lw=lw)
ax1.axhline(0.0, color="k", linestyle="--", lw=lw)
for i in range(0, len(ops)):
    i = len(ops) - i - 1

    plot_ring = 30 - ring1 - 1
    pv_par = all_vzs_poly_par[i][plot_ring]
    pn_par = all_nzs_poly_par[i][plot_ring]
    dpv_par = pv_par.deriv(1)
    dpn_par = pn_par.deriv(1)
    dpg_par = pv_par * dpn_par + pn_par * dpv_par

    keep = np.logical_and(ss[plot_ring]>=xlim1, ss[plot_ring]<=xlim2)
    x = ss[plot_ring][keep]
    y = dpg_par(ss[plot_ring][keep])
    #ax1.plot(x, y, label=i, color=colors[i], linestyle=":")

    # Can at least try and estimate the second derivative of the density by
    # using average values in each region between these knots for each ring.
    avg_nzs = all_nzs[i][:, keep].mean(axis=1)
    rad_dist = dists[:, keep].mean(axis=1)
    nz_grad1 = np.gradient(avg_nzs, rad_dist)
    nz_grad2 = np.gradient(grad1, rad_dist)
    avg_dg_rad = -dperps[i] * nz_grad2[plot_ring]
    #ax1.axhline(avg_dg_rad, linestyle="--", color=colors[i])

    #label = r"$\mathdefault{D}_{\perp}$" + " = {:.1f}".format(dperps[i]) + r" $\mathdefault{m^2/s}$"
    label = r"$\mathdefault{D}_{\perp}$" + " = {:.1f}".format(dperps[i])
    ax1.plot(x, y+avg_dg_rad, color=colors[i], lw=lw, label=label)

    # The further criteria.
    zero_idx = np.abs(x - 39).argmin()
    nz_stag = pn_par(39)
    nz_stag = pn_par(x[zero_idx])
    #vzs_lin = all_vzs_lin[i][plot_ring]
    #vz_slope = vzs_lin.coef[1] * vzs_lin.mapparms()[1]

    vz = all_vzs[i][plot_ring][keep]
    popt, pcov = curve_fit(lin_fit, x, vz)
    vz_slope = popt[0]

    fig2, ax2 = plt.subplots()
    ax2.plot(x, vz)
    ax2.plot(x, lin_fit(x, *popt))
    #ax2.plot(x, vzs_lin(x))
    #ax2.plot(x, vzs_lin.coef[1]*(vzs_lin.mapparms()[0] + vzs_lin.mapparms()[1]*x) + vzs_lin.coef[0])
    fig2.tight_layout()
    fig2.show()

    print("dperp = {:.1f}".format(dperps[i]))
    print("nz = {:.2e}".format(nz_stag))
    print("nz_max = {:.2e}".format(-1/vz_slope * avg_dg_rad))
    print()

ax1.legend(fontsize=12, loc=(0.65, 0.15), framealpha=1.0)
ax1.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax1.set_ylabel(r"$\frac{\mathdefault{d}\mathrm{\Gamma_{||}}^\mathdefault{conv}}{\mathdefault{ds}}$ + $\frac{\mathdefault{d}\mathrm{\Gamma_{\perp}}^\mathdefault{diff}}{\mathdefault{dr}}$", fontsize=fontsize+4)
ax1.set_xlim([xlim1, xlim2])
ax1.grid()
#ax1.set_ylim([-100, 10])
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_yscale("symlog", linthresh=1)

# Gotta set the yticks so the negatives show up.
ytick_labels = [r"-$\mathdefault{10^1}$", r"-$\mathdefault{10^0}$", "0", r"$\mathdefault{10^0}$", r"$\mathdefault{10^1}$", r"$\mathdefault{10^2}$"]
ax1.set_yticks([-10, -1, 0, 1, 10, 100])
ax1.set_yticklabels(ytick_labels)

fig.tight_layout()
fig.show()
