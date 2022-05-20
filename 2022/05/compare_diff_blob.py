# Script to generate some plots that compare a diffusive DIVIMP case with a
# blobby one.
import oedge_plots
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams["font.family"] = "Century Gothic"


print("Loading runs...")
diff_path = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-mrc-shifted-drifts-quick.nc"
blob_path = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-div6.nc"
diff = oedge_plots.OedgePlots(diff_path)
blob = oedge_plots.OedgePlots(blob_path)
absfac = diff.absfac

# A plot for me, the difference in W density on the grid.
print("2D differences plot...")
diff_nz = diff.read_data_2d("DDLIMS", charge="all") * absfac
blob_nz = blob.read_data_2d("DDLIMS", charge="all") * absfac
diff.plot_contour_polygon("KTEBS", own_data=blob_nz-diff_nz, normtype="symlog",
    cbar_label="blob-diff (m-3)", cmap="coolwarm", vmin=-1e16, vmax=1e16)

# Radial profiles of the W density at the OMP.
print("Profiles at OMP plot...")
zomp = float(diff.nc["Z0"][:])
romp = 2.260
rstart = 2.22; rend = 2.37
diff_r, diff_nz = diff.fake_probe(rstart, rend, zomp, zomp, "nz", show_plot=False, plot="R")
blob_r, blob_nz = blob.fake_probe(rstart, rend, zomp, zomp, "nz", show_plot=False, plot="R")
dmask = diff_nz > 0
bmask = blob_nz > 0
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(diff_r[dmask]-romp, diff_nz[dmask]*absfac, label="Diffusive", lw=3, color="tab:red")
ax.plot(blob_r[bmask]-romp, blob_nz[bmask]*absfac, label="Blobby", lw=3, color="tab:purple")
ax.legend(fontsize=14)
ax.tick_params(labelsize=12)
ax.grid()
ax.set_xlim([0, 0.1])
ax.set_xlabel(r"R-$\mathdefault{R_{sep}}$ OMP (m)", fontsize=14)
ax.set_ylabel(r"W Density ($\mathdefault{m^{-3}}$)", fontsize=14)
ax.set_yscale("log")
fig.tight_layout()
fig.show()

# Comparisons along a near-SOL ring.
print("Near-SOL plot...")
ring = 21
rmrso = float(diff.nc["MIDIST"][1][ring])
s, diff_nz = diff.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
s, blob_nz = blob.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(s, diff_nz, lw=3, color="tab:red", label="Diffusive")
ax.plot(s, blob_nz, lw=3, color="tab:purple", label="Blobby")
ax.text(0.03, 0.7, r"R-$\mathdefault{R_{sep}}$"+" OMP = {:.2f} mm".format(rmrso * 1000), transform=ax.transAxes, fontsize=12, bbox=dict(facecolor='white', edgecolor="white", alpha=0.5))
ax.set_xlabel("Distance from inner target (m)", fontsize=14)
ax.set_ylabel(r"W Density ($\mathdefault{m^{-3}}$)", fontsize=14)
ax.legend(fontsize=14)
ax.tick_params(labelsize=12)
ax.set_yscale("log")
ax.grid()
fig.tight_layout()
fig.show()

# Density in the core region overlaid with SXR measurements.
print("Core concentration plot...")
aminor = (2.260 - 1.723)
#aminor = 0.588
r0 = 1.723
core_rings = np.arange(5, 20)
diff_cws = []; blob_cws = []; rhos = []
for ring in core_rings:
    s, ne = diff.along_ring(ring, "KNBS", plot_it=False)
    s, diff_nz = diff.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    s, blob_nz = blob.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    diff_cws.append((diff_nz/ne).mean())
    blob_cws.append((blob_nz/ne).mean())

    # Calculate rho.
    s, r = diff.along_ring(ring, "RS", plot_it=False)
    romp = r.max()
    rho = (romp - r0) / aminor
    rhos.append(rho)
diff_cws = np.array(diff_cws)
blob_cws = np.array(blob_cws)

path = "/Users/zamperini/Documents/d3d_work/files/imp_analysis_167196.npz"
imps = np.load(path)
rho = imps["rho"]
cw = imps["cw"]
cw_avg = np.nanmean(cw, axis=0)
cw_std = np.nanstd(cw, axis=0)

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(rho, cw_avg*100, color="k")
ax.fill_between(rho, (cw_avg-cw_std)*100, (cw_avg+cw_std)*100, color="k", alpha=0.3)
ax.plot(rhos, diff_cws*100, lw=3, color="tab:red", label="Diffusive")
ax.plot(rhos, blob_cws*100, lw=3, color="tab:purple", label="Blobby")
ax.legend(fontsize=14)
ax.set_xlim([0.7, 1])
ax.set_ylim([0, 0.5])
ax.set_xlabel(r"$\mathdefault{\rho}$", fontsize=14)
ax.set_ylabel("W Concentration (%)", fontsize=14)
fig.tight_layout()
fig.show()

# Plot comparisons with the collector probe. Unsure what's going on here.
print("DIVIMP style collector probe plot...")
cp_path = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-div6.collector_probe"
cp = pd.read_csv(cp_path, skiprows=7, delim_whitespace=True)
a2_path = "/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
a2 = pd.read_excel(a2_path)
a2_itf_x = a2["R-Rsep omp D (cm)"].values
a2_otf_x = a2["R-Rsep omp U (cm)"].values
a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values * 1e15 * 1e4  # W/cm2 to W/m2
a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values * 1e15 * 1e4  # W/cm2 to W/m2
exposure_time = 4 * 25
cp_itf_x = cp["ROMP"] * 100  # m to cm
cp_itf_y = cp["IMPFLUX_IDF"] * absfac * exposure_time
cp_otf_x = cp["ROMP"] * 100  # m to cm
cp_otf_y = cp["IMPFLUX_ODF"] * absfac * exposure_time
fig, ax = plt.subplots()
ax.plot(cp_itf_x, cp_itf_y)
ax.plot(cp_otf_x, cp_otf_y)
ax.scatter(a2_itf_x, a2_itf_y, marker="*")
ax.scatter(a2_otf_x, a2_otf_y, marker="*")
ax.set_xlim([6, 15])
ax.set_ylim([1e17, 1e22])
ax.set_yscale("log")
fig.tight_layout()
fig.show()
