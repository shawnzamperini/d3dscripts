# Script to compare some parameters from DIVIMP when swapping out the
# entire wall with SiC.
import oedge_plots
import numpy as np


# The 003 set of runs are the one I initially used for the PVR slides (Dperp = 0.6 m2/s I think).
# The 006 set has been developed using a blobby transport model.
# cpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/pre-github/d3d-allsic-wall-c-003.nc"
# spath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/pre-github/d3d-allsic-wall-si-003.nc"
# gpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/pre-github/d3d-sic-wall-allgrap-003.nc"
cpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-009.nc"
spath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-009.nc"
gpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-009.nc"
cop = oedge_plots.OedgePlots(cpath)
sop = oedge_plots.OedgePlots(spath)
gop = oedge_plots.OedgePlots(gpath)


ne2d = cop.read_data_2d("KNBS")
zeff_sic = np.zeros(ne2d.shape)

dsum1_gph = np.zeros(ne2d.shape)
dsum2_gph = np.zeros(ne2d.shape)
dsum1_sic = np.zeros(ne2d.shape)
dsum2_sic = np.zeros(ne2d.shape)

no_core = False
for charge in range(1, 15):
    print(charge)
    if charge < 7:
        tmp = cop.read_data_2d("DDLIMS", scaling=cop.absfac, charge=charge, no_core=no_core)
        dsum1_sic += tmp * charge
        dsum2_sic += tmp * charge**2


        tmp = gop.read_data_2d("DDLIMS", scaling=gop.absfac, charge=charge, no_core=no_core)
        dsum1_gph += tmp * charge
        dsum2_gph += tmp * charge**2

    tmp = sop.read_data_2d("DDLIMS", scaling=sop.absfac, charge=charge, no_core=no_core)
    dsum1_sic += tmp * charge
    dsum2_sic += tmp * charge**2

# Compute Zeff for graphite walls.
zeff_gph = np.zeros(ne2d.shape)
zeffs1 = dsum1_gph
zeffs2 = 1.0 * ne2d - zeffs1
for i in range(0, len(zeff_gph)):
    if zeffs2[i] > 0:
        zeff_gph[i] = (1.0 * zeffs2[i] + dsum2_gph[i]) / (1.0 * ne2d[i])
    else:
        zeff_gph[i] = dsum2_gph[i] / dsum1_gph[i]

# Compute Zeff for SiC walls.
zeff_sic = np.zeros(ne2d.shape)
zeffs1 = dsum1_sic
zeffs2 = 1.0 * ne2d - zeffs1
for i in range(0, len(zeff_sic)):
    if zeffs2[i] > 0:
        zeff_sic[i] = (1.0 * zeffs2[i] + dsum2_sic[i]) / (1.0 * ne2d[i])
    else:
        zeff_sic[i] = dsum2_sic[i] / dsum1_sic[i]

cbar_label = r"$\mathdefault{Z_{eff}^{SiC}-Z_{eff}^{graphite}}$"
cop.plot_contour_polygon("DDLIMS", own_data=zeff_sic-zeff_gph, cbar_label=cbar_label,
    cmap="coolwarm", lut=21, vmin=-1, vmax=1, normtype="symlog")

# Present as a percent change.
cbar_label = r"$\mathdefault{Z_{eff}^{SiC}-Z_{eff}^{graphite}}$ %"
cop.plot_contour_polygon("DDLIMS", own_data=(zeff_sic-zeff_gph)/zeff_gph * 100, cbar_label=cbar_label,
    cmap="coolwarm", lut=21, vmin=-25, vmax=25, normtype="symlin")