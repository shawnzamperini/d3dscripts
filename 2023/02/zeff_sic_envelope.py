# This is a "back of the envelope" type calculation to compare how Zeff between graphite and SiC could change.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../2022/08")
import EngelhardtModel

# Inputs.
# fc_grap = 0.02
# fo = 0.001
fo = 0.0
ne = 1e19

# Hijack this just to get the interpolation functions.
em = EngelhardtModel.EngelhardtModel()
em.load_mm_model()

# Initial goal is a plot of Te vs. Zeff using the SiC yields.
te = np.geomspace(1, 50, 500)
zeff_grap = np.zeros(te.shape)
zeff_sic = np.zeros(te.shape)
fc_graps = np.zeros(te.shape)
for i in range(0, len(te)):

    # Calculate yields and fractions of C and Si from SiC.
    eimp = 5 * te[i]
    fc_grap = em.Y_D_C(eimp) + em.Y_D_C_ch(eimp)
    fc_sic = (fc_grap * em.Y_D_SiC_C(eimp)) / (em.Y_D_C(eimp) + fc_grap * em.Y_C_C(eimp) - fc_grap * em.Y_C_SiC_C(eimp))
    fsi_sic = (em.Y_D_SiC_Si(eimp) + fc_sic * em.Y_C_SiC_Si(eimp)) / (em.Y_D_SiC_C(eimp) + fc_sic * em.Y_C_SiC_C(eimp))

    # Now calculate densities of each.
    ni_grap = ne / (1 + 6 * fc_grap + 8 * fo)
    no_grap = fo * ni_grap
    nc_grap = fc_grap * ni_grap
    ni_sic = ne / (1 + 6 * fc_sic + 8 * fo + 14 * fsi_sic)
    no_sic = fo * ni_sic
    nc_sic = fc_sic * ni_sic
    nsi_sic = fsi_sic * ni_sic

    # Now Zeff.
    zeff_grap[i] = (6**2 * nc_grap + 8**2 * no_grap + ni_grap) / ne
    zeff_sic[i] = (6**2 * nc_sic + 8**2 * no_sic + ni_sic + 14**2 * nsi_sic) / ne
    fc_graps[i] = fc_grap


fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.plot(te, zeff_grap, label="Graphite", color="tab:red", lw=3)
ax1.plot(te, zeff_sic, label="SiC", color="tab:purple", lw=3)
#ax1.ticklabel_format(style="plain", axis="y")
ax1.legend(fontsize=12)
ax1.set_xlim([0, None])
ax1.set_ylim([1, 5])
ax1.set_yticks(np.arange(1, 5))
#ax1.set_yscale("log")
ax1.grid(which="both", alpha=0.3)
ax1.set_xlabel(r"$\mathdefault{T_e}$ (eV)", fontsize=14)
ax1.set_ylabel(r"$\mathdefault{Z_{eff}}$", fontsize=14)
ax1.tick_params(labelsize=12)


fig.tight_layout()
fig.show()