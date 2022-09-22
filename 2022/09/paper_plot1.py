from importlib import reload
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append("../08/")
import EngelhardtModel

plt.rcParams["font.family"] = "Century Gothic"

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



# Plotting.
fig, ax1 = plt.subplots(figsize=(5,4))
ax1.axvline(em.rsep*100, color="k", linestyle="--")
ax1.axvline(em.riz*100, color="k", linestyle="--")
ax1.axvline(em.rlim*100, color="k", linestyle="--")
#ax1.plot(gph["rs"], zeff_sic, color="tab:purple", lw=3, label="SiC")
#ax1.plot(gph["rs"], zeff_gph, color="tab:red", lw=3, label="Graphite")
ax1.plot(gph["rs"]*100, gph["nz"], color="tab:red", lw=3, label="Graphite")
#ax1.plot(gph["rs"]*100, c_sic["nz"], color="tab:purple", lw=3, label="SiC (C)")
#ax1.plot(gph["rs"]*100, si_sic["nz"]*5, color="tab:purple", lw=3, linestyle="--", label="SiC (Si) x 5")
ax1.set_xlabel("Distance from separatrix (cm)", fontsize=14)
#ax1.set_ylabel("Zeff", fontsize=14)
ax1.set_ylabel(r"$\mathdefault{n_z\ (m^{-3})}$", fontsize=14)
ax1.set_xlim([-0.01*100, em.rwall*100])
ax1.set_ylim([0, 1.1e17])
#ax1.set_ylim([1, None])
ax1.legend(fontsize=14)
fig.tight_layout()
fig.show()
