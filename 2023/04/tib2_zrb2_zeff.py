import oedge_plots
import matplotlib.pyplot as plt
import numpy as np


grap = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-009.nc")
op_ti = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/tib2_zrb2/d3d-tib2-wall-ti-009.nc")
op_ti_b = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/tib2_zrb2/d3d-tib2-wall-b-009.nc")
op_zr = 0
op_zr_b = 0

dsum1 = 0.0
dsum2 = 0.0
dsum1_ti = 0.0
dsum2_ti = 0.0
dsum1_ti_b = 0.0
dsum2_ti_b = 0.0
dsum1_c6 = 0.0
dsum2_c6 = 0.0

# Calculate Zeff using the way DIVIMP does it. I admit I don't fully understand it, but obvi DIVIMP is trustworthy.
for charge in range(1, 23):
    print("Charge: {}".format(charge))

    # Boron.
    if charge < 6:
        nz = np.array(op_ti_b.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
        dsum1_ti_b += nz * charge
        dsum2_ti_b += nz * charge ** 2

    # Carbon.
    if charge < 7:
        nz = np.array(grap.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
        dsum1 += nz * charge
        dsum2 += nz * charge ** 2

        # The CER Zeff calculation is assuming C6+ is the only impurity charge
        # state that exists, so that's what we want to compare to.
        if charge == 6:
            dsum1_c6 += nz * charge
            dsum2_c6 += nz * charge ** 2

    # Titanium.
    if charge < 23:
        nz = np.array(op_ti.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
        dsum1_ti += nz * charge
        dsum2_ti += nz * charge ** 2

    # Zirconium.
    # To do.


# Needed for below.
ne_probe = grap.fake_probe(1.9, 2.5, 0.0, 0.0, "ne")
divimp_r = np.array(ne_probe["r"])
ne = np.array(ne_probe["ne"])
zeff_tib2 = np.zeros(ne.shape)
zeff_gph = np.zeros(ne.shape)

# Mimicking calculation from div.f.
zeffs1 = dsum1
zeffs2 = 1.0 * ne - zeffs1
zeff_gph = np.zeros(zeffs1.shape)
for i in range(0, len(zeff_gph)):
    if zeffs2[i] > 0:
        zeff_gph[i] = (1.0 * zeffs2[i] + dsum2[i]) / (1.0 * ne[i])
    else:
        zeff_gph[i] = dsum2[i] / dsum1[i]

# Again for TiB2.
dsum1_tib2 = dsum1_ti_b + dsum1_ti
dsum2_tib2 = dsum2_ti_b + dsum2_ti
zeffs1 = dsum1_tib2
zeffs2 = 1.0 * ne - zeffs1
zeff_tib2 = np.zeros(len(zeffs1))
for i in range(0, len(zeff_tib2)):
    if zeffs2[i] > 0:
        zeff_tib2[i] = (1.0 * zeffs2[i] + dsum2_tib2[i]) / (1.0 * ne[i])
    else:
        zeff_tib2[i] = dsum2_tib2[i] / dsum1_tib2[i]


divmask = divimp_r > 2.21
fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.axvline(2.26, color="k", linestyle="--")
ax1.plot(divimp_r[divmask], zeff_gph[divmask], color="tab:red", label="Graphite", lw=3)
ax1.plot(divimp_r[divmask], zeff_tib2[divmask], color="tab:purple", label="TiB2", lw=3)
ax1.legend()
ax1.legend(fontsize=14)
ax1.set_xlabel("R (m)", fontsize=16)
ax1.set_ylabel(r"$\mathdefault{Z_{eff}}$", fontsize=16)
fig.tight_layout()
fig.show()

