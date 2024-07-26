# Make plots of Zeff along the separatrix for a hypothetical SiC midplane wall.
import oedge_plots
import numpy as np
import matplotlib.pyplot as plt


cpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-c-002.nc"
spath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-si-002.nc"
gpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-grap-002.nc"
ogpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-othergrap-002.nc"
agpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-002.nc"
cop = oedge_plots.OedgePlots(cpath)
sop = oedge_plots.OedgePlots(spath)
gop = oedge_plots.OedgePlots(gpath)
ogop = oedge_plots.OedgePlots(ogpath)
agop = oedge_plots.OedgePlots(agpath)

# Load density.
sepring = 18
s, ne = cop.along_ring(sepring, "KNBS", plot_it=False)
ne2d = cop.read_data_2d("KNBS")

# Get all the impurity densities along the sep ring.
s, nc = cop.along_ring(sepring, "DDLIMS", charge="all", plot_it=False)
ncs = np.zeros((6, len(nc)))
nss = np.zeros((14, len(nc)))
zeff = np.ones(len(nc))
zeff_grap = np.ones(len(nc))
zeff2d = np.ones(ne2d.shape)
zeff_ag = np.ones(ne2d.shape)
zeff_og = np.ones(ne2d.shape)
zeff_ag1d = np.ones(len(nc))
zeff_og1d = np.ones(len(nc))
for charge in range(1, 15):
    print(charge)
    if charge < 7:
        s, tmp1 = cop.along_ring(sepring, "DDLIMS", charge=charge, plot_it=False)
        ncs[charge-1] = tmp1
        zeff += charge**2 * tmp1 / ne
        s, tmp1 = gop.along_ring(sepring, "DDLIMS", charge=charge, plot_it=False)
        #ncs[charge-1] = tmp1
        zeff_grap += charge**2 * tmp1 / ne

        # Same but for 2D.
        tmp3 = cop.read_data_2d("DDLIMS", scaling=cop.absfac, charge=charge)
        zeff2d += charge**2 * tmp3 / ne2d

        # The other graphite wall portions and an all graphite wall.
        s, tmp7 = ogop.along_ring(sepring, "DDLIMS", charge=charge, plot_it=False)
        zeff_og1d += charge**2 * tmp7 / ne
        tmp5 = ogop.read_data_2d("DDLIMS", scaling=ogop.absfac, charge=charge)
        zeff_og += charge**2 * tmp5 / ne2d

        s, tmp8 = agop.along_ring(sepring, "DDLIMS", charge=charge, plot_it=False)
        zeff_ag1d += charge**2 * tmp8 / ne
        tmp6 = agop.read_data_2d("DDLIMS", scaling=agop.absfac, charge=charge)
        zeff_ag += charge**2 * tmp6 / ne2d

    s, tmp2 = sop.along_ring(sepring, "DDLIMS", charge=charge, plot_it=False)
    #nss[charge-1] = tmp2
    zeff += charge**2 * tmp2 / ne

    # Same but for 2D.
    tmp4 = sop.read_data_2d("DDLIMS", scaling=sop.absfac, charge=charge)
    zeff2d += charge**2 * tmp4 / ne2d

fig, ax = plt.subplots()
ax.plot(s, zeff+zeff_og1d-1, label="SiC vertical plate")
ax.plot(s, zeff_ag1d, label="All graphite")
ax.set_xlabel("Distance from inner target (m)", fontsize=14)
ax.set_ylabel("Zeff @ separatrix", fontsize=14)
ax.legend()
fig.tight_layout()
fig.show()

# Contribution due only to wall portion.
cop.plot_contour_polygon("DDLIMS", own_data=zeff2d, cbar_label="Zeff (SiC plate only)",
    cmap="nipy_spectral", lut=21, vmin=1, vmax=1.05, normtype="log")

# Full Zeff (graphite + SiC vertical plate).
cop.plot_contour_polygon("DDLIMS", own_data=zeff2d+zeff_og-1, cbar_label="Zeff (total)",
    cmap="nipy_spectral", lut=21, vmin=1, vmax=5, normtype="log")

# Full Zeff (all graphite wall).
cop.plot_contour_polygon("DDLIMS", own_data=zeff_ag, cbar_label="Zeff (all graphite)",
    cmap="nipy_spectral", lut=21, vmin=1, vmax=5, normtype="log")

# Zeff change with vertical SiC plate.
cbar_label = r"$\mathdefault{Z_{eff}^{SiC}-Z_{eff}^{graphite}}$"
cop.plot_contour_polygon("DDLIMS", own_data=(zeff2d+zeff_og-1)-zeff_ag, cbar_label=cbar_label,
    cmap="coolwarm", lut=21, vmin=-0.5, vmax=0.5, normtype="symlin")
