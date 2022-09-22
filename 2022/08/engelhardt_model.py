from importlib import reload
import EngelhardtModel
import matplotlib.pyplot as plt
import numpy as np


# For testing, easy to just do this here.
reload(EngelhardtModel)

# Inputs.
riz            = 0.07
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


em = EngelhardtModel.EngelhardtModel()
em.set_geometry(riz=riz, rlim=rlim, rwall=rwall)
em.set_ne(nesep=nesep, lambda_ne=lambda_ne)
em.set_te_ti(tesep=tesep, lambda_te=lambda_te, timult=timult)
em.load_mm_model()
out1 = em.run_mm_model()

# Run separately for graphite and SiC (both atoms).
dperp = [1, 1, 5, 10]
conns = [100, 100, 50, 10]
vels = [500, 500, 500, 500]
c_sic = em.run_engelhardt(mat="c_sic", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)
si_sic = em.run_engelhardt(mat="si_sic", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)
gph = em.run_engelhardt(mat="c", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)

# Calculate the Zeff profiles for graphite and SiC.
zeff_sic = em.calc_zeff_prof_sic(c_sic["rs"], c_sic["nz"], si_sic["nz"])
zeff_gph = em.calc_zeff_prof_gph(gph["rs"], gph["nz"])


# Plotting.
fig, ax1 = plt.subplots(figsize=(5,4))
ax1.axvline(em.rsep, color="k", linestyle="--")
ax1.axvline(em.riz, color="k", linestyle="--")
ax1.axvline(em.rlim, color="k", linestyle="--")
#ax1.plot(gph["rs"], zeff_sic, color="tab:purple", lw=3, label="SiC")
#ax1.plot(gph["rs"], zeff_gph, color="tab:red", lw=3, label="Graphite")
ax1.plot(gph["rs"], c_sic["nz"], color="tab:purple", lw=3, label="SiC")
ax1.plot(gph["rs"], gph["nz"], color="tab:red", lw=3, label="Graphite")
ax1.set_xlabel("R-Rsep (m)", fontsize=14)
#ax1.set_ylabel("Zeff", fontsize=14)
ax1.set_ylabel("nz (m-3)", fontsize=14)
ax1.set_xlim([-0.01, em.rwall])
#ax1.set_ylim([1, None])
ax1.legend()
fig.tight_layout()
fig.show()

# Do a scan in lambda_ne and riz to see how it affects Zeff at the separatrix.
n = 20
lamb_nes = np.linspace(0.02, 0.13, n)
rizs = np.linspace(0.05, 0.08, n)
lamb_zeffs_sic = []
lamb_zeffs_gph = []
lamb_nzs_sic_si = []
lamb_nzs_sic_c = []
lamb_nzs_gph = []
riz_zeffs_sic = []
riz_zeffs_gph = []
riz_nzs_sic_si = []
riz_nzs_sic_c = []
riz_nzs_gph = []
for i in range(0, n):

    # Lambda_ne first.
    em.set_geometry(riz=riz, rlim=rlim, rwall=rwall)
    em.set_ne(nesep=nesep, lambda_ne=lamb_nes[i])
    em.set_te_ti(tesep=tesep, lambda_te=lambda_te, timult=timult)
    em.run_mm_model()
    c_sic = em.run_engelhardt(mat="c_sic", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)
    si_sic = em.run_engelhardt(mat="si_sic", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)
    gph = em.run_engelhardt(mat="c", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)
    zeff_sic = em.calc_zeff_prof_sic(c_sic["rs"], c_sic["nz"], si_sic["nz"])
    zeff_gph = em.calc_zeff_prof_gph(gph["rs"], gph["nz"])
    lamb_zeffs_sic.append(zeff_sic[0])
    lamb_zeffs_gph.append(zeff_gph[0])
    lamb_nzs_sic_c.append(c_sic["nz"][0])
    lamb_nzs_sic_si.append(si_sic["nz"][0])
    lamb_nzs_gph.append(gph["nz"][0])

    # Then rizs.
    em.set_geometry(riz=rizs[i], rlim=rlim, rwall=rwall)
    em.set_ne(nesep=nesep, lambda_ne=lambda_ne)
    em.set_te_ti(tesep=tesep, lambda_te=lambda_te, timult=timult)
    em.run_mm_model()
    c_sic = em.run_engelhardt(mat="c_sic", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)
    si_sic = em.run_engelhardt(mat="si_sic", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)
    gph = em.run_engelhardt(mat="c", dperp=dperp, conns=conns, match_peter=match_peter, mode=mode, vels=vels)
    zeff_sic = em.calc_zeff_prof_sic(c_sic["rs"], c_sic["nz"], si_sic["nz"])
    zeff_gph = em.calc_zeff_prof_gph(gph["rs"], gph["nz"])
    riz_zeffs_sic.append(zeff_sic[0])
    riz_zeffs_gph.append(zeff_gph[0])
    riz_nzs_sic_c.append(c_sic["nz"][0])
    riz_nzs_sic_si.append(si_sic["nz"][0])
    riz_nzs_gph.append(gph["nz"][0])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4), sharey=True)
ax1.plot(lamb_nes, lamb_nzs_sic_c, label="SiC (C)", color="tab:purple", lw=3)
ax1.plot(lamb_nes, lamb_nzs_sic_si, label="SiC (Si)", color="tab:purple", linestyle="--", lw=3)
ax1.plot(lamb_nes, lamb_nzs_gph, label="Graphite", color="tab:red", lw=3)
ax2.plot(rizs, riz_nzs_sic_c, label="SiC (C)", color="tab:purple", lw=3)
ax2.plot(rizs, riz_nzs_sic_si, label="SiC (Si)", color="tab:purple", linestyle="--", lw=3)
ax2.plot(rizs, riz_nzs_gph, label="Graphite", color="tab:red", lw=3)
ax1.set_ylabel(r"$\mathdefault{n_z\ @\ r_{iz}\ (m^{-3})}$", fontsize=14)
ax1.set_xlabel(r"$\mathdefault{\lambda_{ne}\ (m)}$", fontsize=14)
ax2.set_xlabel(r"$\mathdefault{r_{iz}\ (m)}$", fontsize=14)
ax1.grid(alpha=0.3)
ax2.grid(alpha=0.3)
ax2.legend(fontsize=14)
fig.tight_layout()
fig.show()

# A comparison to Peter's example.
