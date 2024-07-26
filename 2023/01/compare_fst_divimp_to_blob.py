import oedge_plots
import matplotlib.pyplot as plt
import LimPlots
import numpy as np
from scipy.interpolate import griddata
import pickle


print("Loading DIVIMP runs...")
# Drifts OFF. Dperp = 0.3. Similar to used in FST paper but this uses a differnet background (d3d-167196-bg-shifted-ring-entry-10).
# fst_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mrc-shifted-nodrift-2-copy.nc"
fst_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-diff-001-predep.nc"
# Drifts ON (60%). Dperp = 0.3. bkg 010.
#fst_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mrc-shifted-drifts-3a.nc"
fst = oedge_plots.OedgePlots(fst_nc_path)
# Drifts ON. bkg 010.
# blob_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-004.nc"
# Drifts ON. bkg 13. tcorr = 5 us.
# blob_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-010f-predep.nc"
blob_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-fluc-002.nc"
blob = oedge_plots.OedgePlots(blob_nc_path)
absfac = fst.absfac

ring = 91
s, fst_nz = fst.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
s, blob_nz = blob.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
fst_nz = fst_nz / fst_nz.max()
blob_nz = blob_nz / blob_nz.max()

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(s, fst_nz, label="FST Paper")
ax.plot(s, blob_nz, label="Blobby")
ax.legend()
ax.set_xlabel("S (m)")
ax.set_ylabel("W Density (normalized)")
fig.tight_layout()
fig.show()

# Similar to how I did in tomas_sxr_results.py, I will compare the radial profiles from DIVIMP vs the one from 3DLIM
# to see how they line up in magnitude. I will compare the normal vs. the blobby DIVIMP ones above to a 3DLIM run where
# the necessary absfac has been determined based off what gives the best match to the profiles in magnitude.
fst_probe = fst.fake_probe(2.21, 2.36, -0.188, -0.188, "nz", charge="all")
blob_probe = blob.fake_probe(2.21, 2.36, -0.188, -0.188, "nz", charge="all")

# Load the gfile and map the R, Z CER points to psin.
print("Creating psin interpolation...")
gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]

def get_lim_data(lim_path, lim_absfac):
    lp = LimPlots.LimPlots(lim_path)
    lpdata = lp.plot_par_rad("nz", 21, charge="all")

    # Let's not choose right in the middle (Y = 0) since flows are zero there, and it actually leads to distracting
    # valleys in the radial profile due to no impurities there.
    mididx = np.argmin(np.abs(lpdata["X"][:, 0])) - 5

    rorigin = 2.282
    rad_locs = rorigin - lpdata["Y"][0]

    # This value is chosen in lim_a2_again.py. It is that needed to match the experimental deposition.
    lim_nz = lpdata["Z"][mididx].data * lim_absfac
    lim_rzs = zip(rad_locs, np.full(len(rad_locs), -0.188))
    lim_psins = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(lim_rzs))
    lim_rhos = np.sqrt(lim_psins)
    wall_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), (rorigin + 0.08, -0.188))
    wall_rho = np.sqrt(wall_psin)
    return {"psin":lim_psins, "nz":lim_nz}

lim_fst = get_lim_data("/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240_44-noprobe.nc", 1e17)
#lim_blob = get_lim_data("/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240_44-blob-comp-001a.nc", 1e17)
lim_blob = get_lim_data("/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-011-noprobe.nc", 1e17)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True, sharey=True)
ax2.plot(blob_probe["psin"], np.array(blob_probe["nz"])*absfac, label="Blobby")
ax2.plot(lim_blob["psin"], lim_blob["nz"], label="3DLIM")
ax1.plot(fst_probe["psin"], np.array(fst_probe["nz"])*absfac, label="FST")
ax1.plot(lim_fst["psin"], lim_fst["nz"], label="3DLIM")
ax1.set_yscale("log")
ax1.legend()
ax1.grid(alpha=0.3)
ax2.grid(alpha=0.3)
ax1.set_xlabel("Psin")
ax2.set_xlabel("Psin")
ax1.set_ylabel("W Density (m-3)")
fig.tight_layout()
fig.show()


# div_mask = np.array(blob_probe["psin"]) < 1.14
div_mask = np.array(blob_probe["psin"]) < 9999
div_psin = np.array(blob_probe["psin"])[div_mask]
div_blob_nz = np.array(blob_probe["nz"])[div_mask]*absfac
div_fst_nz = np.array(fst_probe["nz"])[div_mask]*absfac
fig, ax = plt.subplots(figsize=(5, 4))
ax.axvline(1.0, color="k")
ax.plot(div_psin, div_blob_nz, label="Blobby", color="tab:red", lw=3)
ax.plot(div_psin, div_fst_nz, label="Diffusive", color="tab:purple", lw=3)
ax.plot(lim_blob["psin"], lim_blob["nz"], label="3DLIM", color="tab:cyan", lw=3)
ax.set_yscale("log")
ax.legend(fontsize=14)
ax.grid(alpha=0.3)
ax.set_xlabel(r"$\mathdefault{\psi_N}$", fontsize=14)
ax.set_ylabel(r"W Density $\mathdefault{(m^{-3})}$", fontsize=14)
fig.tight_layout()
fig.show()