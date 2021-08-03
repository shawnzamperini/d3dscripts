from scipy.interpolate import Rbf
import numpy as np
import oedge_plots
import matplotlib.pyplot as plt
from matplotlib import ticker
from numpy.polynomial import Polynomial


dperp = 0.3
path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc"
op = oedge_plots.OedgePlots(path)

# Ring number - 1 for indexing.
ring1 = 16 - 1; ring2 = 45 - 1
num_charges = 30

sep_r = op.nc.variables["RS"][:][ring1]
sep_z = op.nc.variables["ZS"][:][ring1]

def get_rmrs(r, z):
    dists = np.sqrt((sep_r-r)**2 + (sep_z-z)**2)
    return dists[dists.argmin()]

# Create an Rbf interpolation of thetag and distance from the separatrix.
thetags = np.zeros((ring2-ring1, len(sep_r)))
rmrss = np.zeros(thetags.shape)
vzs = np.zeros((num_charges, thetags.shape[0], thetags.shape[1]))
for ring in range(0, ring2-ring1):
    print("Indexing ring {}".format(ring+ring1))
    ring_thetag = op.nc.variables["THETAG"][:][ring+ring1]
    ring_rmrs = np.zeros(len(ring_thetag))
    for knot in range(0, len(ring_thetag)):
        r = op.nc.variables["RS"][:][ring+ring1, knot]
        z = op.nc.variables["ZS"][:][ring+ring1, knot]
        ring_rmrs[knot] = get_rmrs(r, z)
    thetags[ring] = ring_thetag
    rmrss[ring] = ring_rmrs

    # Append the parallel velocities.
    for charge in range(0, num_charges):

        # Additional +1 bc we subtracted 1 at the top of the script but this function
        # will want the actual ring number.
        s, vz = op.along_ring(ring+ring1+1, "VELAVG", charge=charge, plot_it=False, remove_zeros=False)
        vzs[charge, ring] = vz

# Net density.
nzs = op.nc.variables["DDLIMS"][:][1:-1].sum(axis=0)[ring1:ring2].data
nzs_charge = op.nc.variables["DDLIMS"][:][1:-1][:num_charges, ring1:ring2].data

print("  Calculating weighted averages...")
num_rings = thetags.shape[0]; num_knots = thetags.shape[1]
vzs_weighted = np.zeros((num_rings, num_knots))
nzs_summed   = np.zeros((num_rings, num_knots))
nzs_poly_par = []; vzs_poly_par = []; vzs_lin = []
for i in range(0, num_rings):
    for j in range(0, num_knots):
        nzs_summed[i, j] = nzs_charge[:, i, j].sum()
        if nzs_summed[i, j] != 0.0:
            vzs_weighted[i, j] = np.average(vzs[:, i, j], weights=nzs_charge[:, i, j])

# Remove all knots where thetag = 0.
ss = op.nc.variables["KSS"][:][ring1:ring2].data
keep_cols = thetags[0] != 0
thetags = thetags[:, keep_cols]
rmrss = rmrss[:, keep_cols]
nzs = nzs[:, keep_cols]
ss = ss[:, keep_cols]
vzs_weighted = vzs_weighted[:, keep_cols]

rbf = Rbf(thetags, rmrss, nzs, epsilon=0.01)

nzs_int = rbf(thetags, rmrss)

#nz_min = min(nzs.min(), nzs_int.min())
nz_min = 0.0001
nz_max = max(nzs.max(), nzs_int.max())
levels = np.linspace(nz_min, nz_max, 10)

nzs_masked = np.ma.masked_where(nzs<=0, nzs)
nzs_int_masked = np.ma.masked_where(nzs_int<=0, nzs_int)
#nzs_masked_ds1 = np.gradient(nzs_masked, thetags, axis=1)
#nzs_masked_dr1 = np.gradient(nzs_masked, rmrss, axis=0)

buffer = 0.05  # in meters, radially
s_buffer = 5  # in meters, parallel
num_points = 31
rad_nzs_d1s = np.zeros(thetags.shape)
rad_nzs_d2s = np.zeros(thetags.shape)
par_nzs_d1s = np.zeros(thetags.shape)
par_vzs_d1s = np.zeros(thetags.shape)
par_gz_d1s = np.zeros(thetags.shape)
rad_gz_d1s = np.zeros(thetags.shape)
for i in range(0, thetags.shape[0]):

    # Fit an Rbf to the parallel velocity and density.
    rbf_nzpar = Rbf(ss[i], nzs[i])
    rbf_vzpar = Rbf(ss[i], vzs_weighted[i])

    # Fit to the nonzero values.
    nonzero = vzs_weighted[i] != 0.0
    try:
        lin_vzpar = Polynomial(1).fit(ss[i][nonzero], vzs_weighted[i][nonzero], 1)
    except:
        print("Ring {} (i={}) vzs are all zeros. Exiting.".format(i+ring1+1, i))

    for j in range(0, thetags.shape[1]):
        tg = np.full(num_points, thetags[i, j])
        rm = np.linspace(max(rmrss[i, j]-buffer, 0), rmrss[i, j]+buffer, num_points)
        rad_nz = rbf(tg, rm)
        rad_nz_d1 = np.gradient(rad_nz, rm)
        rad_nz_d2 = np.gradient(rad_nz_d1, rm)

        # Return the middle point.
        rad_nzs_d1s[i, j] = rad_nz_d1[int(num_points/2)]
        rad_nzs_d2s[i, j] = rad_nz_d2[int(num_points/2)]

        # Parallel derivatives.
        s = np.linspace(max(ss[i, j]-s_buffer, 0), ss[i, j]+s_buffer, num_points)
        nzpar_int = rbf_nzpar(s)
        #vzpar_int = rbf_vzpar(s)
        vzpar_int = lin_vzpar(s)
        nzpar_d1 = np.gradient(nzpar_int, s)
        vzpar_d1 = np.gradient(vzpar_int, s)

        # Return the middle point again.
        par_nzs_d1s[i, j] = nzpar_d1[int(num_points/2)]
        par_vzs_d1s[i, j] = vzpar_d1[int(num_points/2)]
        par_gz_d1s[i, j] = nzpar_int[int(num_points/2)] * vzpar_d1[int(num_points/2)] + vzpar_int[int(num_points/2)] * nzpar_d1[int(num_points/2)]
        rad_gz_d1s[i, j] = - dperp * rad_nz_d2[int(num_points/2)]

# Contour plot of the density.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

#cont1 = ax1.contourf(thetags, rmrss, nzs_masked, cmap="magma", locator=ticker.LogLocator(), levels=levels)
cont1 = ax1.contourf(thetags, rmrss, nzs_masked, cmap="magma")
cont2 = ax2.contourf(thetags, rmrss, par_gz_d1s+rad_gz_d1s, cmap="coolwarm", levels=np.linspace(-10, 10, 11), extend="both")
fig.colorbar(cont1, ax=ax1)
fig.colorbar(cont2, ax=ax2)

fig.tight_layout()
fig.show()

# Tag on a plot for PFMC of the velocity along a ring with the density for all
# charge states.
plot_ring = 18 - 16 - 1
vz_y = vzs_weighted[plot_ring] / np.max(np.abs(vzs_weighted[plot_ring]).max())
nz_y = nzs[plot_ring] / np.max(np.abs(nzs[plot_ring]).max())

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 5

fig, ax1 = plt.subplots()
ax1.axhline(0.0, color="k", linestyle="--", lw=lw)
ax1.plot(ss[plot_ring], vz_y, label=r"$\mathdefault{v_z}$", lw=lw, color=colors[1])
ax1.plot(ss[plot_ring], nz_y, label=r"$\mathdefault{n_z}$", lw=lw, color=colors[3])
#ax1.legend(fontsize=fontsize)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax1.set_ylabel("Normalized values", fontsize=fontsize)
ax1.grid()
ax1.set_xlim([20, 60])
ax1.set_ylim([-0.5, None])
ax1.tick_params(axis='both', which='major', labelsize=12)
fig.tight_layout()
fig.show()
