# For 167196 convert psin's from the LPs to ring numbers.
import pandas as pd
import oedge_plots
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt


# Load needed data.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d-167196-modE-shelf-bg-shifted-smooth-2.nc"
op = oedge_plots.OedgePlots(ncpath)
xlpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/data_for_167196.xlsx"
df = pd.read_excel(xlpath, sheet_name="lp_data")
df = df[['psin', 'r', 'z', 'pname', 'te', 'ne (cm-3)', 'jsat (A/cm2)',
       'ne (m-3)', 'jsat (A/m2)', 'shot', 'target']]

# Do not want to worry about inner target data since we only had one probe.
mask2 = df["target"] == "outer"
df = df[mask2]

# Represent the LP data as interpolations of a smoothed signal.
df = df.sort_values("psin")
df["te_smooth"] = savgol_filter(df["te"], 15, 2)
df["ne_smooth"] = savgol_filter(df["ne (m-3)"], 15, 2)
f_te = interp1d(df["psin"], df["te_smooth"], bounds_error=False, fill_value=(df["te_smooth"].iloc[0], df["te_smooth"].iloc[-1]))
f_ne = interp1d(df["psin"], df["ne_smooth"], bounds_error=False, fill_value=(df["ne_smooth"].iloc[0], df["ne_smooth"].iloc[-1]))

# Will compare to grid psins to find closest ring index.
psitar = op.nc["PSITAR"][0]
rings = np.arange(1, len(psitar)+1)

# Need to exclude core rings. Handled separately.
mask = rings >= 19
psitar = psitar[mask]
rings = rings[mask]

# Also create an interpolation that maps a psin value to a ring.
f_ring = interp1d(psitar, rings)

lp_psins = df["psin"]
lp_rings = [round(float(f_ring(psin))) for psin in lp_psins]
#for psin in lp_psins:
#    ring = round(f_ring(psin))
#    idx = np.argmin(np.abs(psin-psitar))
#    ring = rings[idx]
#    lp_rings.append(ring)

# We want to give a ring (all of them from 19 to 190) and have it give us the
# respective ne and Te. This is facilitated by an "alternate" ring that goes
# from 0-179, which is in order of psin. This makes an interpolation possible.
# LP: te --> psin --> alt_ring(psin) --> ring(alt_ring)
# LP: ring --> alt_ring(ring) --> psin(alt_ring) --> te(psin)
sort_idx = np.argsort(psitar)
ring_sort = rings[sort_idx]
psitar_sort = psitar[sort_idx]
alt_rings = np.arange(0, len(psitar_sort))
ring_to_alt = dict(zip(ring_sort, alt_rings))
f_alt_psin = interp1d(alt_rings, psitar_sort)
ring_vals = {"ring":[], "psin":[], "te":[], "ne":[]}
for ring in range(19, 191):
    alt_ring = ring_to_alt[ring]
    psin = float(f_alt_psin(alt_ring))
    te = float(f_te(psin))
    ne = float(f_ne(psin))
    ring_vals["ring"].append(ring)
    ring_vals["psin"].append(psin)
    ring_vals["te"].append(te)
    ring_vals["ne"].append(ne)
ring_df = pd.DataFrame(ring_vals)

# Put into dataframe, then sort and average for each ring.
#df["lp_rings"] = lp_rings
#ring_df = df.groupby("lp_rings").mean().sort_index()
ring_df.to_excel("tmp2.xlsx")

# Plotting to see if it makes sense.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5), sharex=True)
ax1.scatter(df["psin"], df["te"], alpha=0.6, color="k")
ax1.plot(ring_df["psin"], ring_df["te"], color="r", lw=2)
ax2.scatter(df["psin"], df["ne (m-3)"], alpha=0.6, color="k")
ax2.plot(ring_df["psin"], ring_df["ne"], color="r", lw=2)
ax1.set_ylabel("Te eV)")
ax2.set_ylabel("ne (m-3)")
ax1.set_xlim([0.97, 1.1])
fig.supxlabel("Psin")
fig.tight_layout()
fig.show()
