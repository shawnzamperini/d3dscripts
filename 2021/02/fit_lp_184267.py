import get_lp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Inputs for the fitting function.
lp_xl_path = "/mnt/c/Users/Shawn/Documents/d3d_work/lp_data/lps_184267_271-272_v1.xlsx"
lp_xl_sheet = "Data Fixed Psin Sorted"

# Fit the jsat data.
jsat_dict = get_lp.fit_conv_gauss(lp_xl_path, lp_xl_sheet,
  lp_xl_ydata="jsat fixed (A/cm2)", gauss_range=[1.0, 1.04],
  ylabel="jsat (A/cm2)")
jsat_psin = jsat_dict["psin_fit"]
jsat_fit  = jsat_dict["y_fit"]

# Fit the Te data.
te_dict = get_lp.fit_conv_gauss(lp_xl_path, lp_xl_sheet,
  lp_xl_ydata="Te fixed (eV)", gauss_range=[1.0, 1.04],
  ylabel="Te (eV)")
te_psin = te_dict["psin_fit"]
te_fit  = te_dict["y_fit"]

# Plot to see what the end result really is.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
ax1.plot(jsat_dict["psin"], jsat_dict["y"], '.', color="k", ms=3)
ax1.plot(jsat_psin, jsat_fit, color="k", lw=4)
ax1.plot(jsat_psin, jsat_fit, color="tab:red", lw=3)
ax1.set_ylabel("jsat (A/cm2)", fontsize=16)
ax1.set_xlabel("psin", fontsize=16)
ax2.plot(te_dict["psin"], te_dict["y"], '.', color="k", ms=3)
ax2.plot(te_psin, te_fit, color="k", lw=4)
ax2.plot(te_psin, te_fit, color="tab:red", lw=3)
ax2.set_ylabel("Te (eV)", fontsize=16)
ax2.set_xlabel("psin", fontsize=16)
ax1.grid()
ax2.grid()
fig.tight_layout()
fig.show()

# Save the output in the format needed for the DIVIMP input file.
empties = np.zeros(len(jsat_psin))
df = pd.DataFrame({"Psin":jsat_psin, "Te":te_fit, "Ti":te_fit,
  "jsat":jsat_fit*1000, "empty1":empties, "empty2":empties, "empty3":empties})
df.to_excel("184267_271-272_lp_outer_fit.xlsx")
