# This script is to compare radial profiles from OSM to the data from Thomson
# and a reciprocating Langmuir probe.
from oedge_plots import OedgePlots
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter


# A note of Mach directions. For our set of USN DIVIMP runs, S=0 is at the
# outer target, so a positive Mach number is towards the inner target. For the
# RCP data, to line up with this convention, we need to multiply the RCP Mach
# numbers by -1 for 184527.

# Constants and inputs.
shot = 184527

if shot == 184527:
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/184527/d3d-184527-bkg-004.nc"
    ts_path = "/Users/zamperini/Documents/d3d_work/divimp_files/184527/omfit_184527_ts_v4.xlsx"
    core_shift = 0.0
    div_shift = 0.0
    mach_mult = -1
elif shot == 184267:
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/184267/d3d-184267-bkg-003.nc"
    ts_path = "/Users/zamperini/Documents/d3d_work/divimp_files/184267/omfit_184267_ts.xlsx"
    core_shift = 0.02
    div_shift = 0.02
    mach_mult = 1
rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"


def calc_er(x, vp):
    """
    Can eventually write a function to smooth or fit Vp to get Er from.
    """

    vp_sg = savgol_filter(vp, 11, 2)
    er = np.gradient(vp_sg, x)
    return er


# Load everything in.
op = OedgePlots(ncpath)
ts = pd.read_excel(ts_path)
rcp = pd.ExcelFile(rcp_path)
mp1 = rcp.parse("MP{}_1".format(shot))
mp2 = rcp.parse("MP{}_2".format(shot))
xp1 = rcp.parse("XP{}_1".format(shot))
xp2 = rcp.parse("XP{}_2".format(shot))

# Print out the psin of each ring in the SOL.
psins = op.nc.variables["PSIFL"][:,0]
print("{:4} {:.4}".format("Ring", "Psin"))
for n, p in enumerate(psins):
    if p == 0.0:
        continue
    print("{:4} {:.4f}".format(n+1, p))

# Extract all the data we wish to plot, put into convienent dictionaries.
# o = osm, t = ts and r = rcp.
o = {}; t = {}; r = {}

# OSM data along the RCP plunges and TS locations first.
mp_rmin = 2.25;  mp_rmax = 2.36;  mp_zmin = -0.188; mp_zmax = -0.188  # MiMES
xp_rmin = 1.485; xp_rmax = 1.485; xp_zmin = -1.25;  xp_zmax = -1.01   # XP
cr_rmin = 1.94;  cr_rmax = 1.94;  cr_zmin =  0.58;  cr_zmax =  0.82   # Core
dv_rmin = 1.485; dv_rmax = 1.485; dv_zmin = -1.25;  dv_zmax = -1.01   # Divertor, same as XP
sa_rmin = 1.485; sa_rmax = 1.485; sa_zmin =  0.93;  sa_zmax =  1.18   # SAS
print("Loading fake MiMES data...")
o["m_psin"], o["m_te"] = op.fake_probe(mp_rmin, mp_rmax, mp_zmin, mp_zmax, "Te",   plot="psin", show_plot=False, verbal=False)
o["m_psin"], o["m_ne"] = op.fake_probe(mp_rmin, mp_rmax, mp_zmin, mp_zmax, "ne",   plot="psin", show_plot=False, verbal=False)
o["m_psin"], o["m_m"]  = op.fake_probe(mp_rmin, mp_rmax, mp_zmin, mp_zmax, "Mach", plot="psin", show_plot=False, verbal=False)
o["m_psin"], o["m_er"] = op.fake_probe(mp_rmin, mp_rmax, mp_zmin, mp_zmax, "Erad", plot="psin", show_plot=False, verbal=False)
#o["mp_psin"], o["mp_ep"] = op.fake_probe(mp_rmin, mp_rmax, mp_zmin, mp_zmax, "Epol", plot="psin", show_plot=False, verbal=False)
print("Loading fake XP data...")
o["x_psin"], o["x_te"] = op.fake_probe(xp_rmin, xp_rmax, xp_zmin, xp_zmax, "Te",   plot="psin", show_plot=False, verbal=False)
o["x_psin"], o["x_ne"] = op.fake_probe(xp_rmin, xp_rmax, xp_zmin, xp_zmax, "ne",   plot="psin", show_plot=False, verbal=False)
o["x_psin"], o["x_m"]  = op.fake_probe(xp_rmin, xp_rmax, xp_zmin, xp_zmax, "Mach", plot="psin", show_plot=False, verbal=False)
o["x_psin"], o["x_er"] = op.fake_probe(xp_rmin, xp_rmax, xp_zmin, xp_zmax, "Erad", plot="psin", show_plot=False, verbal=False)
#o["xp_psin"], o["xp_ep"] = op.fake_probe(xp_rmin, xp_rmax, xp_zmin, xp_zmax, "Epol", plot="psin", show_plot=False, verbal=False)
print("Loading fake core TS data...")
o["c_psin"], o["c_te"] = op.fake_probe(cr_rmin, cr_rmax, cr_zmin, cr_zmax, "Te",   plot="psin", show_plot=False, verbal=False)
o["c_psin"], o["c_ne"] = op.fake_probe(cr_rmin, cr_rmax, cr_zmin, cr_zmax, "ne",   plot="psin", show_plot=False, verbal=False)
print("Loading fake divertor TS data...")
o["d_psin"], o["d_te"] = o["x_psin"], o["x_te"]
o["d_psin"], o["d_ne"] = o["x_psin"], o["x_ne"]
print("Loading fake SAS TS data...")
o["s_psin"], o["s_te"] = op.fake_probe(sa_rmin, sa_rmax, sa_zmin, sa_zmax, "Te",   plot="psin", show_plot=False, verbal=False)
o["s_psin"], o["s_ne"] = op.fake_probe(sa_rmin, sa_rmax, sa_zmin, sa_zmax, "ne",   plot="psin", show_plot=False, verbal=False)

# Thomson data next.
core = ts[ts["System"] == "core"]
div = ts[ts["System"] == "divertor"]
sas = ts[ts["System"] == "divertor_sas"]
t["c_psin"] = core["Psin"].values + core_shift
t["c_te"]   = core["Te (eV)"].values
t["c_ne"]   = core["ne (m-3)"].values
t["d_psin"] = div["Psin"].values + div_shift
t["d_te"]   = div["Te (eV)"].values
t["d_ne"]   = div["ne (m-3)"].values
t["s_psin"] = sas["Psin"].values
t["s_te"]   = sas["Te (eV)"].values
t["s_ne"]   = sas["ne (m-3)"].values

# RCP data last.
r["mp1_psin"] = mp1["Psin"]
r["mp1_te"] = mp1["Te (eV)"]
r["mp1_ne"] = mp1["ne (1e18 m-3)"] * 1e18
r["mp1_m"]  = mp1["Mach"] * mach_mult
r["mp1_vp"] = (mp1["Vf2 (V)"] + mp1["Vf3 (V)"]) / 2.0
r["mp1_x"] = mp1["R (cm)"] / 100
r["mp2_psin"] = mp2["Psin"]
r["mp2_te"] = mp2["Te (eV)"]
r["mp2_ne"] = mp2["ne (1e18 m-3)"] * 1e18
r["mp2_m"] = mp2["Mach"] * mach_mult
r["mp2_vp"] = (mp2["Vf2 (V)"] + mp2["Vf3 (V)"]) / 2.0
r["mp2_x"] = mp2["R (cm)"] / 100
r["xp1_psin"] = xp1["Psin"]
r["xp1_te"] = xp1["Te (eV)"]
r["xp1_ne"] = xp1["ne (1e18 m-3)"] * 1e18
r["xp1_m"] = xp1["Mach"] * mach_mult
r["xp1_x"] = xp1["Z (cm)"] / 100
r["xp2_psin"] = xp2["Psin"]
r["xp2_te"] = xp2["Te (eV)"]
r["xp2_ne"] = xp2["ne (1e18 m-3)"] * 1e18
r["xp2_m"] = xp2["Mach"] * mach_mult
r["xp2_x"] = xp2["Z (cm)"] / 100

# Er from the RCPs needs to be calculated.
r["mp1_er"] = calc_er(r["mp1_x"], r["mp1_vp"])
r["mp2_er"] = calc_er(r["mp2_x"], r["mp2_vp"])


def plot_comparison(ax, osm_x, osm_y, exp_x, exp_y, exp_x2=None, exp_y2=None,
    drop_vals=[], label=None, ylabel=None, xlims=None, ylims=None, logy=False):
    """
    Given an Axes object, plot OSM data to compare to the given experimental
    data at the same location.
    """

    # Drop nans and zeros.
    keep = ~np.isnan(osm_x)
    x = osm_x[keep]
    y = osm_y[keep]
    keep = ~np.isnan(osm_y)
    x = osm_x[keep]
    y = osm_y[keep]
    for val in drop_vals:
        keep = y != val
        x = x[keep]
        y = y[keep]

    ax.scatter(exp_x, exp_y, s=10, marker="^", alpha=0.7, color="tab:red", edgecolor="k")
    # Can have a second set of experimental data if more than one plunge.
    if type(exp_x2) != type(None):
        ax.scatter(exp_x2, exp_y2, s=10, marker="^", alpha=0.7, color="tab:red", edgecolor="k")
    ax.plot(x, y, color="tab:red", label=label)
    ax.axvline(1.0, color="k", linestyle="--", lw=1)
    ax.set_ylabel(ylabel, fontsize=10)
    #ax.legend(fontsize=8)
    if type(xlims) != type(None):
        ax.set_xlim(xlims)
    if type(ylims) != type(None):
        ax.set_ylim(ylims)
    if logy:
        ax.set_yscale("log")

# Now for a large grid of plots showing comparisons with experimental data. I
# Expect somewhere near 15 plots, but each is pretty basic.
fig, axs = plt.subplots(4, 5, figsize=(12, 8))
axs = axs.flatten()

for i in range(0, len(axs)):
    axs[i].tick_params(labelsize=8)

axs[0].set_title("MiMES RCP", fontsize=10)
axs[1].set_title("XP RCP", fontsize=10)
axs[2].set_title("Core TS", fontsize=10)
axs[3].set_title("Divertor TS", fontsize=10)
axs[4].set_title("SAS TS", fontsize=10)

# First row Te
te_xlims = [0.98, 1.15]
te_ylims = [0, 100]
plot_comparison(axs[0], o["m_psin"], o["m_te"], r["mp1_psin"], r["mp1_te"], r["mp2_psin"], r["mp2_te"], drop_vals=[0, 1], label="MRCP", ylabel="Te (eV)", xlims=te_xlims, ylims=te_ylims)
plot_comparison(axs[1], o["x_psin"], o["x_te"], r["xp1_psin"], r["xp1_te"], r["xp2_psin"], r["xp2_te"], drop_vals=[0, 1], label="XRCP", xlims=te_xlims, ylims=te_ylims)
plot_comparison(axs[2], o["c_psin"], o["c_te"], t["c_psin"],   t["c_te"],  drop_vals=[0, 1], label="Core", xlims=te_xlims, ylims=te_ylims)
plot_comparison(axs[3], o["d_psin"], o["d_te"], t["d_psin"],   t["d_te"],  drop_vals=[0, 1], label="Divertor", xlims=te_xlims, ylims=te_ylims)
plot_comparison(axs[4], o["s_psin"], o["s_te"], t["s_psin"],   t["s_te"],  drop_vals=[0, 1], label="SAS", xlims=te_xlims, ylims=te_ylims)

# Second row ne
ne_xlims = te_xlims
ne_ylims = [0, 2e19]
plot_comparison(axs[5], o["m_psin"], o["m_ne"], r["mp1_psin"], r["mp1_ne"], r["mp2_psin"], r["mp2_ne"], drop_vals=[0, o["m_ne"][36]], label="MRCP", ylabel="ne (m-3)", ylims=ne_ylims)
plot_comparison(axs[6], o["x_psin"], o["x_ne"], r["xp1_psin"], r["xp1_ne"], r["xp2_psin"], r["xp2_ne"], drop_vals=[0, o["m_ne"][36]], label="XRCP", ylims=ne_ylims)
plot_comparison(axs[7], o["c_psin"], o["c_ne"], t["c_psin"],   t["c_ne"],  drop_vals=[0, o["m_ne"][36]], label="Core", xlims=ne_xlims, ylims=ne_ylims)
plot_comparison(axs[8], o["d_psin"], o["d_ne"], t["d_psin"],   t["d_ne"],  drop_vals=[0, o["m_ne"][36]], label="Divertor", xlims=ne_xlims, ylims=ne_ylims)
plot_comparison(axs[9], o["s_psin"], o["s_ne"], t["s_psin"],   t["s_ne"],  drop_vals=[0, o["m_ne"][36]], label="SAS", xlims=ne_xlims, ylims=ne_ylims)

# Third row Mach
m_xlims = te_xlims
m_ylims = [-0.5, 0.5]
plot_comparison(axs[10], o["m_psin"], o["m_m"],  r["mp1_psin"], r["mp1_m"],  r["mp2_psin"], r["mp2_m"],  drop_vals=[0], label="MRCP", ylabel="Mach", xlims=m_xlims, ylims=m_ylims)
plot_comparison(axs[11], o["x_psin"], o["x_m"],  r["xp1_psin"], r["xp1_m"],  r["xp2_psin"], r["xp2_m"],  drop_vals=[0], label="XRCP", xlims=m_xlims, ylims=m_ylims)

# Fourth row Er
er_xlims = te_xlims
er_ylims = None
plot_comparison(axs[15], o["m_psin"], o["m_er"],  r["mp1_psin"], r["mp1_er"],  r["mp2_psin"], r["mp2_er"],  drop_vals=[0], label="MRCP", ylabel="Er (V/m)", xlims=er_xlims, ylims=er_ylims, logy=False)

fig.tight_layout()
fig.show()
