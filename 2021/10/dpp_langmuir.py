# Script to make a pretty plot of the Langmuir probe data for 184527.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


plt.rcParams['font.family'] = 'Century Gothic'
#plt.rc('axes', unicode_minus=False)

xl_path = "/Users/zamperini/Documents/d3d_work/184527/lp_184527.xlsx"
df = pd.read_excel(xl_path, sheet_name="For Python")

psin_out     = df["Psin"][df["Target"] == "Outer"]
te_out       = df["Te (eV)"][df["Target"] == "Outer"]
psin_out_fit = df["Fit Psin"][:99]
te_out_fit   = df["Fit Te (eV)"][:99]

fig, ax = plt.subplots(figsize=(5, 4))

ax.scatter(psin_out, te_out, s=30, color="tab:red", edgecolor="k", label="Measured")
ax.plot(psin_out_fit, te_out_fit, color="k", lw=2, label="Fit")
ax.axvline(1, color="k", linestyle="--", lw=2)
ax.set_xlabel("Psin", fontsize=16)
ax.set_ylabel("Te (eV)", fontsize=16)
ax.grid()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xticks(np.arange(0.99, 1.03, 0.01))

fig.tight_layout()
fig.show()
