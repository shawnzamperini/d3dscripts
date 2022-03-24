import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

xl_path = "/Users/zamperini/My Drive/Research/Documents/2022/01/lp_167196.xlsx"
df = pd.read_excel(xl_path, skiprows=1)

# Data is stored in the first 6 columns.
data = df[df.columns[:6]].dropna()

# Additionally drop anything with negative densities.
data = data[data["ne (1e18 m-3)"] > 0]

r = data["R (cm)"].values
te = data["Te (eV)"].values
ne = data["ne (1e18 m-3)"].values

# Create bins and assign values to each.
bins = np.linspace(r.min(), r.max(), 20)
bin_idxs = np.digitize(r, bins)

# Compute binned averages.
r_binned = []
te_binned = []
ne_binned = []
for bin_idx in np.unique(bin_idxs):
    bin_r = r[bin_idxs==bin_idx].mean()
    bin_te = te[bin_idxs==bin_idx].mean()
    bin_ne = ne[bin_idxs==bin_idx].mean()
    r_binned.append(bin_r)
    te_binned.append(bin_te)
    ne_binned.append(bin_ne)
r_binned = np.array(r_binned)
te_binned = np.array(te_binned)
ne_binned = np.array(ne_binned)

# Print out for copy/paste to Excel file.
print("R (cm)  Te (eV)  ne (1e18 m-3)")
for i in range(0, len(r_binned)):
    print("{:.2f}\t{:.2f}\t{:.2f}".format(r_binned[::-1][i], te_binned[::-1][i], ne_binned[::-1][i]))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

ax1.scatter(r, te, color="r", s=20, marker="^", edgecolors=None, alpha=0.3)
ax1.scatter(r_binned, te_binned, color="r", s=40, marker="^", edgecolors="k")
ax1.set_yscale("log")
ax1.set_xlabel("R (cm)")
ax1.set_ylabel("Te (eV)")

ax2.scatter(r, ne)
ax2.scatter(r_binned, ne_binned)
ax2.set_yscale("log")
ax2.set_xlabel("R (cm)")
ax2.set_ylabel("ne (1e18 m-3)")


fig.tight_layout()
fig.show()
