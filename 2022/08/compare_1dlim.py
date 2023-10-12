# Script to compare the profiles from a 1DLIM run to that of the simple
# Engelhardt model.
from importlib import reload
import EngelhardtModel
import matplotlib.pyplot as plt
import numpy as np
import LimPlots

plt.rcParams["font.family"] = "Century Gothic"

# For testing, easy to just do this here.
reload(EngelhardtModel)

# Inputs.
riz_c          = 0.07; riz_si = 0.095
rlim           = 0.10
rwall          = 0.12
nesep          = 1e19
tesep          = 100
lambda_ne      = 0.05
lambda_te      = 0.05
timult         = 3.0
mat            = "c_sic"
match_peter    = False
mode           = "diffusive"

# The 1DLIM paths.
# 001: SiC, C - diff
# 002: SiC, Si - diff
# 003: SiC, Si - conv
# 004: SiC, C - conv
ncpath_c = "/Users/zamperini/Documents/d3d_work/lim_runs/tor_lim_testing/engelhardt-1dlim-001.nc"
ncpath_si = "/Users/zamperini/Documents/d3d_work/lim_runs/tor_lim_testing/engelhardt-1dlim-002.nc"


em = EngelhardtModel.EngelhardtModel()


# Run separately for graphite and SiC (both atoms).
dperp = [1, 1, 1, 1]
conns = [100, 100, 100, 100]
vels = [500, 500, 500, 500]

# SiC, C
em.set_geometry(riz=riz_c, rlim=rlim, rwall=rwall)
em.set_ne(nesep=nesep, lambda_ne=lambda_ne)
em.set_te_ti(tesep=tesep, lambda_te=lambda_te, timult=timult)
em.load_mm_model()
out1 = em.run_mm_model()
c_sic = em.run_engelhardt(mat="c_sic", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)

# SiC, Si
em.set_geometry(riz=riz_si, rlim=rlim, rwall=rwall)
em.set_ne(nesep=nesep, lambda_ne=lambda_ne)
em.set_te_ti(tesep=tesep, lambda_te=lambda_te, timult=timult)
em.load_mm_model()
out1 = em.run_mm_model()
si_sic = em.run_engelhardt(mat="si_sic", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)

# Graphite
em.set_geometry(riz=riz_c, rlim=rlim, rwall=rwall)
em.set_ne(nesep=nesep, lambda_ne=lambda_ne)
em.set_te_ti(tesep=tesep, lambda_te=lambda_te, timult=timult)
em.load_mm_model()
out1 = em.run_mm_model()
gph = em.run_engelhardt(mat="c", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)

# Calculate the Zeff profiles for graphite and SiC.
zeff_sic = em.calc_zeff_prof_sic(c_sic["rs"], c_sic["nz"], si_sic["nz"])
zeff_gph = em.calc_zeff_prof_gph(gph["rs"], gph["nz"])

# Load 1DLIM run. Change R coordinates to be R-Rsep.
lp_c = LimPlots.LimPlots(ncpath_c)
nz_dict_c = lp_c.plot_par_rad("nz", 20, charge="all", showplot=False)
lp_nc = nz_dict_c["Z"].data[0]
lp_r = 0.10 - nz_dict_c["Y"][0]

# Normalize the 1DLIM data to 1 at the limiter location, and then multiply
# by what we know the density should be there from the mixed-material model above.
lim_idx = np.where(np.abs(lp_r-0.10)==np.abs(lp_r-0.10).min())
lp_fact = lp_nc[lim_idx[0][0]]
lp_nc = lp_nc / lp_fact
fact_sic_c = out1["fc_sic"] * em.fne(0.10)
lp_nc = lp_nc * fact_sic_c

# Function from LimPlots on shortening the data to what matters.
def trim(Z):
    xouts = lp_c.nc['XOUTS'][:].data
    youts = lp_c.nc['YOUTS'][:].data
    xwids = lp_c.nc["XWIDS"][:].data
    mask = xwids != 0
    xouts = xouts[mask]
    xwids = xwids[mask]
    xkeep_min = np.nonzero(xouts)[0].min()
    xkeep_max = np.nonzero(xouts)[0].max()
    ykeep_min = np.nonzero(youts)[0].min()
    ykeep_max = np.nonzero(youts)[0].max()
    xouts = xouts[xkeep_min:xkeep_max]
    youts = youts[ykeep_min:ykeep_max]
    xwids = xwids[xkeep_min:xkeep_max]
    Z = Z[ykeep_min:ykeep_max, xkeep_min:xkeep_max]
    yabsorb1a = float(lp_c.nc["yabsorb1a"][:].data)
    yabsorb2a = float(lp_c.nc["yabsorb2a"][:].data)
    ykeep = np.where(np.logical_and(youts>=yabsorb2a, youts<=yabsorb1a))[0]
    youts = youts[ykeep]
    Z = Z[ykeep, :]
    return Z

# Calculate Zeff. Need to pull DDLIM3 so we can do charge-state resolved calcs.
# Scale the value to what it is at the limiter location like before. We need to
# multiply by the absfac calculation to be consistent with the above scaling
# procedure since it is baked into plot_par_rad, messy, I know.
absfac = float(lp_c.nc["ABSFAC"][:].data)
ddlim3 = lp_c.nc["DDLIM3"][20, :, :, :] * absfac / lp_fact
ddlim3_all = trim(ddlim3[1:].sum(axis=0))
ne = lp_c.nc["CRNBS"][:].data
ne = trim(ne)[0]
for charge in range(1, 7):

    # Get the (unscaled) density for just this charge, trim it to just the
    # simulation volume.
    nc_charge = ddlim3[charge]
    nc_charge = trim(nc_charge)[0]

    # Represent the density as what fraction of the total density it is, and
    # then multiply by our scaled density lp_nc.
    charge_frac = nc_charge / ddlim3_all
    nc_charge = nc_charge * lp_nc

    if charge == 1:
        zeff = np.ones(nc_charge.shape)
        zeff_c = np.ones(nc_charge.shape)
    zeff_c += charge**2 * nc_charge / ne
    zeff += charge**2 * nc_charge / ne

# Repeat the above, but with SiC, Si
lp_si = LimPlots.LimPlots(ncpath_si)
nz_dict_si = lp_si.plot_par_rad("nz", 20, charge="all", showplot=False)
lp_nsi = nz_dict_si["Z"].data[0]
#lp_r = 0.10 - nz_dict_c["Y"][0]

#lim_idx = np.where(np.abs(lp_r-0.10)==np.abs(lp_r-0.10).min())
lp_fact = lp_nsi[lim_idx[0][0]]
lp_nsi = lp_nsi / lp_fact
fact_sic_si = out1["fsi_sic"] * em.fne(0.10)
lp_nsi = lp_nsi * fact_sic_si

absfac = float(lp_si.nc["ABSFAC"][:].data)
ddlim3 = lp_si.nc["DDLIM3"][20, :, :, :] * absfac / lp_fact
ddlim3_all = trim(ddlim3[1:].sum(axis=0))
#ne = lp_si.nc["CRNBS"][:].data
#ne = trim(ne)[0]
for charge in range(1, 15):

    # Get the (unscaled) density for just this charge, trim it to just the
    # simulation volume.
    nsi_charge = ddlim3[charge]
    nsi_charge = trim(nsi_charge)[0]

    # Represent the density as what fraction of the total density it is, and
    # then multiply by our scaled density lp_nc.
    charge_frac = nsi_charge / ddlim3_all
    nsi_charge = nsi_charge * lp_nsi

    if charge == 1:

        # We do zeros here since C already includes it. We will sorta add these
        # together later, so don't want to double count.
        zeff_si = np.zeros(nsi_charge.shape)
    zeff_si += charge**2 * nsi_charge / ne
    zeff += charge**2 * nsi_charge / ne

# Do again for just graphite (can reuse the SiC-C case above, just different
# scaling factor is all).
lp_gph = LimPlots.LimPlots(ncpath_c)
nz_dict_gph = lp_gph.plot_par_rad("nz", 20, charge="all", showplot=False)
lp_ngph = nz_dict_gph["Z"].data[0]
#lp_r = 0.10 - nz_dict_c["Y"][0]

#lim_idx = np.where(np.abs(lp_r-0.10)==np.abs(lp_r-0.10).min())
lp_fact = lp_ngph[lim_idx[0][0]]
lp_ngph = lp_ngph / lp_fact
fact_gph = out1["fc_gph"] * em.fne(0.10)
lp_ngph = lp_ngph * fact_gph

absfac = float(lp_gph.nc["ABSFAC"][:].data)
ddlim3 = lp_gph.nc["DDLIM3"][20, :, :, :] * absfac / lp_fact
ddlim3_all = trim(ddlim3[1:].sum(axis=0))
#ne = lp_si.nc["CRNBS"][:].data
#ne = trim(ne)[0]
for charge in range(1, 7):

    # Get the (unscaled) density for just this charge, trim it to just the
    # simulation volume.
    ngph_charge = ddlim3[charge]
    ngph_charge = trim(ngph_charge)[0]

    # Represent the density as what fraction of the total density it is, and
    # then multiply by our scaled density lp_nc.
    charge_frac = ngph_charge / ddlim3_all
    ngph_charge = ngph_charge * lp_ngph

    if charge == 1:

        # We do zeros here since C already includes it. We will sorta add these
        # together later, so don't want to double count.
        zeff_gph = np.ones(ngph_charge.shape)
    zeff_gph += charge**2 * ngph_charge / ne

em_rc = c_sic["rs"]
em_nc = c_sic["nz"]
em_rsi = si_sic["rs"]
em_nsi = si_sic["nz"]

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4), sharex=True)
#fig1, ax1 = plt.subplots(figsize=(5,4))
ax1.axvline(rlim*100, color="k", linestyle="--")
#ax1.axvline(riz_c*100, color="k", linestyle="--")
ax1.plot(lp_r*100, lp_nc, label="SiC", lw=3, color="tab:purple")
#ax1.plot(em_rc*100, em_nc, label="Engelhardt", lw=3, color="tab:red")
ax1.plot(lp_r*100, lp_ngph, label="Graphite", lw=3, color="r")
ax2.axvline(rlim*100, color="k", linestyle="--")
#ax2.axvline(riz_si*100, color="k", linestyle="--")
ax2.plot(lp_r*100, lp_nsi, label="1DLIM", lw=3, color="tab:purple")
#ax2.plot(em_rsi*100, em_nsi, label="Engelhardt", lw=3, color="tab:red")
ax1.legend(fontsize=14, loc="upper left")
ax1.set_xlim([0, 0.12*100])
fig1.supxlabel("Distance from separatrix (cm)", fontsize=16)
ax1.set_title(r"$\mathdefault{Carbon\ density\ (m^{-3})}$", fontsize=16)
ax2.set_title(r"$\mathdefault{Silicon\ density\ (m^{-3})}$", fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax1.set_ylim([0, 1.5e17])
ax2.set_ylim(0, 5e15)
fig1.tight_layout()
fig1.show()

fig2, ax2 = plt.subplots(figsize=(5,4))

# I think this should be zero-based. 4/20/22.
# ax2.fill_between(lp_r*100, 1, zeff_c, color="tab:pink", label="SiC (C)")
# ax2.fill_between(lp_r*100, zeff_c, zeff_c+zeff_si, color="tab:cyan", label="SiC (Si)")
# ax2.plot(lp_r*100, zeff, color="k", lw=3)
# ax2.plot(lp_r*100, zeff_gph, color="r", lw=3, label="Graphite")
ax2.fill_between(lp_r*100, 1-1, zeff_c-1, color="tab:pink", label="SiC (C)")
ax2.fill_between(lp_r*100, zeff_c-1, zeff_c+zeff_si-1, color="tab:cyan", label="SiC (Si)")
ax2.plot(lp_r*100, zeff-1, color="k", lw=3)
ax2.plot(lp_r*100, zeff_gph-1, color="r", lw=3, label="Graphite")
ax2.set_xlabel("Distance from separatrix (cm)", fontsize=16)
ax2.set_ylabel(r"$\mathdefault{Contribution\ to\ Z_{eff}}$", fontsize=16)
ax2.set_xlim([0, 0.12*100])
ax2.set_ylim([1-1, None])
ax2.legend(fontsize=14)
ax2.tick_params(axis='both', which='major', labelsize=12)
fig2.tight_layout()
fig2.show()
