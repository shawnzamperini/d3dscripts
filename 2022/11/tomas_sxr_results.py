#
import pickle
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.interpolate import griddata, interp1d
from gadata import gadata
from scipy.signal import savgol_filter
import oedge_plots
import LimPlots

sys.path.append("/usr/local/mdsplus/python")
import MDSplus

print("Loading ImpRad data...")
path = "/Users/zamperini/Documents/d3d_work/files/190422_w_imprad.pickle"
with open(path, "rb") as f:
    sxr = pickle.load(f)

t = sxr.t.data
rho = sxr.rho.data
wdens = sxr["W dens."].data

# Subset of times to plot radial profiles for.
times = [3000]
wdens_times = []
wdens_err_times = []
for time in times:
    minidx = np.argmin(np.abs(t - time))
    wdens_val = np.array([wdens[minidx, i].n for i in range(0, len(rho))])
    wdens_err = np.array([wdens[minidx, i].s for i in range(0, len(rho))])
    wdens_times.append(wdens_val)
    wdens_err_times.append(wdens_err)

# Load CER data, we want to compare W density to C.
print("Loading CER data...")
conn = MDSplus.Connection("atlas.gat.com")
shot = 190422
tmin = 2800
tmax = 5000
cer = {}
tchords = ["T{}".format(i) for i in range(1, 49)]
vchords = ["V{}".format(i) for i in range(1, 33)]
chords = np.append(tchords, vchords)
for chord in chords:
    print("Chord: {}".format(chord))
    signal_r = "CERAR{}".format(chord)
    signal_z = "CERAZ{}".format(chord)
    signal_zeff = "CERAZEFF{}".format(chord)
    signal_nz = "CERANZ{}".format(chord)
    gaobj_r = gadata(signal_r, shot, connection=conn)
    gaobj_z = gadata(signal_z, shot, connection=conn)
    gaobj_zeff = gadata(signal_zeff, shot, connection=conn)
    gaobj_nz = gadata(signal_nz, shot, connection=conn)
    chord_r = float(gaobj_r.zdata[0])
    chord_z = float(gaobj_z.zdata[0])

    # Return average Zeff value.
    mask = np.logical_and(gaobj_zeff.xdata > tmin, gaobj_zeff.xdata < tmax)
    chord_zeff = float(gaobj_zeff.zdata[mask].mean())
    chord_nz = float(gaobj_nz.zdata[mask].mean())
    chord_zeff_std = float(gaobj_zeff.zdata[mask].std())
    chord_nz_std = float(gaobj_nz.zdata[mask].std())

    # If zero, then no data for this chord.
    if chord_zeff != 0 and not np.isnan(chord_zeff):
        # Time dependent data.
        time = np.array(gaobj_zeff.xdata[mask])
        nz_t = np.array(gaobj_nz.zdata[mask])
        zeff_t = np.array(gaobj_zeff.zdata[mask])

        cer[chord] = {"r": chord_r, "z": chord_z, "zeff": chord_zeff,
                      "nz": chord_nz, "zeff_err": chord_zeff_std, "nz_err": chord_nz_std,
                      "time": time, "nz_t": nz_t, "zeff_t": zeff_t}

# Load the gfile and map the R, Z CER points to psin.
print("Creating psin interpolation...")
gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190423/190423_3000.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]
# tck = bisplrep(Rs, Zs, psin, s=5)

cer_rzs = {}
for chord in cer.keys():
    cer_rzs[chord] = (cer[chord]["r"], cer[chord]["z"])
cer_psins = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(cer_rzs.values()))
count = 0
for chord in cer.keys():
    cer[chord]["psin"] = cer_psins[count]
    count += 1

# Grab data for a chord near rho = 0.2 (this would be CER chord T1, rho ~ 0.18).
# For SXR this is time index 12. Call this the "comp" location.
wdens_comp = np.array([wdens[i, 12].n for i in range(0, len(t))])
wdens_comp_err = np.array([wdens[i, 12].s for i in range(0, len(t))])
f_wdens_comp = interp1d(t, wdens_comp)
cer_comp_t = cer["T1"]["time"]
cer_comp_nz = savgol_filter(cer["T1"]["nz_t"], 21, 1)
wdens_comp_int = f_wdens_comp(cer_comp_t)
avg_nw_nc = (wdens_comp_int / cer_comp_nz).mean()
print("Average nW/nC = {:.2e}".format(avg_nw_nc))

# Now we make a BIG assumption that the ratio of nW/nC stays constant in the
# core, and then make a new profile of W density from the CER data.
nw_cer_est_rho = []
nw_cer_est = []
for chord in cer.keys():
    if chord[0] == "T":
        cer[chord]["nw_est"] = cer[chord]["nz"] * avg_nw_nc
        nw_cer_est_rho.append(np.sqrt(cer[chord]["psin"]))
        nw_cer_est.append(cer[chord]["nz"] * avg_nw_nc)
sort_idx = np.argsort(nw_cer_est_rho)
nw_cer_est_rho = np.array(nw_cer_est_rho)[sort_idx]
nw_cer_est = np.array(nw_cer_est)[sort_idx]

# Load a DIVIMP case so that we can overlay that on our master plot as well.
# tungsten-002: ring-entry mach flow, 0.3 m2/s + blobby @ 1.5 kHz
# tungsten-003: ring-entry mach flow, 0.3 m2/s + blobby @ 3.0 kHz
# tungsten-004: ring-entry mach flow, 0.03 m2/s + blobby @ 3.0 kHz
# tunsgten-005: ring-entry mach flow, 0.00 m2/s + blobby from 167196-blobby-004 @ 2 kHz (Gaussian blob vr at 800 m/s, 400 m/s width)
# tungsten-006: same as 005, but drifts ON.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-tungsten-006.nc"
op = oedge_plots.OedgePlots(ncpath)
div_data = op.fake_probe(2.20, 2.30, 0.0, 0.0, "nz")
div_rho = np.sqrt(div_data["psin"])
nonan = ~np.isnan(div_rho)
div_rho = div_rho[nonan]
absfac = op.absfac
div_nw = np.array(div_data["nz"])[nonan] * absfac

# An additional comparison run without the blobby contribution.
# tungsten-001: ring-entry mach flow, 0.3 m2/s
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-tungsten-001.nc"
op = oedge_plots.OedgePlots(ncpath)
div_data_diff = op.fake_probe(2.20, 2.30, 0.0, 0.0, "nz")
div_rho_diff = np.sqrt(div_data_diff["psin"])
nonan_diff = ~np.isnan(div_rho_diff)
div_rho_diff = div_rho_diff[nonan_diff]
absfac_diff = op.absfac
div_nw_diff = np.array(div_data_diff["nz"])[nonan_diff] * absfac_diff

# Load a 3DLIM run, pull some radial profiles from different locations to 
# capture 3D effects.
lim_path = "/Users/zamperini/Documents/d3d_work/lim_runs/190423/190423-noprobe-002.nc"
lp = LimPlots.LimPlots(lim_path)
lpdata = lp.plot_par_rad("nz", 21, charge="all")

# Let's not choose right in the middle (Y = 0) since flows are zero there, and it actually leads to distracting
# valleys in the radial profile due to no impurities there.
mididx = np.argmin(np.abs(lpdata["X"][:, 0])) - 15

rorigin = 2.305
rad_locs = rorigin - lpdata["Y"][0]
lim_absfac = 2e16  # This value is chosen in lim_mcp01w.py. It is that needed to match the experimental deposition.
lim_nz = lpdata["Z"][mididx].data * lim_absfac
lim_rzs = zip(rad_locs, np.full(len(rad_locs), -0.188))
lim_psins = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(lim_rzs))
lim_rhos = np.sqrt(lim_psins)
wall_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), (rorigin + 0.08, -0.188))
wall_rho = np.sqrt(wall_psin)

# Create a mask for the SXR data so we don't plot the garbage data.
sxr_mask = rho < 0.2
cer_mask = np.logical_and(nw_cer_est_rho >= 0.18, nw_cer_est_rho < 1.015)
div_mask = div_rho > 1.0

fig, ax1 = plt.subplots(figsize=(7, 4))
ax1.axvline(1.0, linestyle="--", color="k")
for i in range(0, len(times)):
    ax1.fill_between(rho[sxr_mask], wdens_times[i][sxr_mask] - wdens_err_times[i][sxr_mask],
                     wdens_times[i][sxr_mask] + wdens_err_times[i][sxr_mask], alpha=0.35, color="tab:red")
    ax1.plot(rho[sxr_mask], wdens_times[i][sxr_mask], label="SXR", color="tab:red", lw=3)
ax1.plot(nw_cer_est_rho[cer_mask], nw_cer_est[cer_mask], label="CER Estimate", color="tab:purple", lw=3)
ax1.plot(div_rho[div_mask][:-1], div_nw[div_mask][:-1], label="DIVIMP - Blobby", color="tab:pink", lw=3)
ax1.plot(div_rho_diff[div_mask][:-1], div_nw_diff[div_mask][:-1], label="DIVIMP - Diffusive", color="tab:pink", lw=3,
         linestyle="--")
ax1.plot(lim_rhos, lim_nz, label="3DLIM", color="tab:green", lw=3)
#ax1.axvline(wall_rho, color="k", lw=2)
ax1.axvline(1.185, color="k", lw=2)
ax1.legend(loc="lower left")
ax1.grid()
ax1.set_yscale("log")
ax1.set_xlim([-0.05, 1.25])
ax1.set_ylim(5e11, 1e16)
# ax1.set_ylim(0, 1.5e1`5)
ax1.set_xlabel("Rho")
ax1.set_ylabel("nW (m-3)")
fig.tight_layout()
fig.show()

fig, (ax2, ax3) = plt.subplots(1, 2, figsize=(9, 4))
ax2.fill_between(t, wdens_comp - wdens_comp_err, wdens_comp + wdens_comp_err, alpha=0.25, color="tab:red")
ax2.plot(t, wdens_comp, color="tab:red")
ax2.set_ylabel("SXR W Density (m-3)", color="tab:red")
ax2.set_xlabel("Time (ms)")
ax2.tick_params(axis="y", labelcolor="tab:red")

ax22 = ax2.twinx()
ax22.plot(cer_comp_t, cer_comp_nz, color="tab:purple")
ax22.set_ylabel("CER C6+ Density (m-3)", color="tab:purple")
ax22.tick_params(axis="y", labelcolor="tab:purple")

ax3.plot(cer_comp_t, wdens_comp_int / cer_comp_nz)
ax3.set_ylabel("nW/nC")
ax3.set_xlabel("Time (ms)")
ax3.set_title("nW/nC @ rho = 0.2")

fig.tight_layout()
fig.show()
