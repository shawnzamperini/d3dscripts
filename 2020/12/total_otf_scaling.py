import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# Which scaling to plot. Parameters used in each option:
# 1: RMRSOMP
# 2: PSOL, TRI
# 3: PSOL, TRI, RMRSOMP
# 4: PSOL
# 5: PSOL RMRSOMP
# 6: PSOL, LAMBDA, RMRSOMP
# 7: PSOL, LCONN
# 8: PSOL, LCONN, RMRSOMP
# 9: LCONN
# 10: PSOL, Q95
plot_scaling = 7
labels = False
plot_none = False
skip_a37 = True
normalize = False

# Path to the Excel file with the data I compiled for Forward Bt, H-mode probes.
xl_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2020/12/w_cp_trends.xlsx"
df = pd.read_excel(xl_path, sheet_name="for_import")

# Pull out data, and ignore the last 3 rows since those are for A13, 14 and 15
# and there isn't data for them, but they would be included if we did.
rows = 13
#print("Using only first 10 probes. Fix when you get the rest of the L data.")
#rows = 9  # Temporary until I get the rest of the connection length data.
probes = df["Probe"].values[:rows]
otf_tot = df["OTF Total/Shot"].values[:rows]
psol = df["PSOL"].values[:rows]
tri = df["Triangularity"].values[:rows]
rmr = df["OTF RMRSOMP Tip (cm)"].values[:rows]
lne = df["Lambda ne (cm)"].values[:rows]
lcon = df["Avg L (m)"].values[:rows]
q95 = df["q95"].values[:rows]

# Normalize values if asked to.
if normalize:
    psol = psol / psol.max()
    tri  = tri  / tri.max()
    rmr  = rmr  / rmr.max()
    lne  = lne  / lne.max()
    lcon = lcon / lcon.max()
    q95  = q95  / q95.max()

# Mask for shelf/floor.
floor_mask = df["Floor/Shelf"].values[:rows] == "Floor"
shelf_mask = df["Floor/Shelf"].values[:rows] == "Shelf"

# Drop A37 since it's such a wonky case.
if skip_a37:
    keep_idx = np.where(probes != "A37")[0]
    probes = probes[keep_idx]
    otf_tot = otf_tot[keep_idx]
    psol = psol[keep_idx]
    tri = tri[keep_idx]
    rmr = rmr[keep_idx]
    lne = lne[keep_idx]
    floor_mask = floor_mask[keep_idx]
    shelf_mask = shelf_mask[keep_idx]

# This is data for probes w/o RBS data that we can plot as holes. Assumed it's the
# rest of the dataframe.
probes_none = df["Probe"].values[rows:]
psol_none = df["PSOL"].values[rows:]
tri_none = df["Triangularity"].values[rows:]
rmr_none = df["OTF RMRSOMP Tip (cm)"].values[rows:]
floor_none = df["Floor/Shelf"].values[rows:] == "Floor"
shelf_none = df["Floor/Shelf"].values[rows:] == "Shelf"

def scaling1(x1, a, b):
    return a * x1**b

def scaling2(X1X2, a, b, c):
    return a * X1X2[0]**b * X1X2[1]**c

def scaling3(X1X2X3, a, b, c, d):
    return a * X1X2X3[0]**b * X1X2X3[1]**c * X1X2X3[2]**d

# Do a simple scaling law.
popt1, pcov1 = curve_fit(scaling1, rmr, otf_tot)
popt2, pcov2 = curve_fit(scaling2, (psol, tri), otf_tot)
popt3, pcov3 = curve_fit(scaling3, (psol, tri, rmr), otf_tot)
popt4, pcov4 = curve_fit(scaling1, psol, otf_tot)
popt5, pcov5 = curve_fit(scaling2, (psol, rmr), otf_tot)
popt6, pcov6 = curve_fit(scaling3, (psol, lne, rmr), otf_tot)
popt7, pcov7 = curve_fit(scaling2, (psol, lcon), otf_tot)
popt8, pcov8 = curve_fit(scaling3, (psol, lcon, rmr), otf_tot)
popt9, pcov9 = curve_fit(scaling1, lcon, otf_tot)
popt10, pcov10 = curve_fit(scaling2, (psol, q95), otf_tot)

# Pull powers for each.
const1    = popt1[0]
rmr_pow1  = popt1[1]
const2    = popt2[0]
psol_pow2 = popt2[1]
tri_pow2  = popt2[2]
const3    = popt3[0]
psol_pow3 = popt3[1]
tri_pow3  = popt3[2]
rmr_pow3  = popt3[3]
const4    = popt4[0]
psol_pow4 = popt4[1]
const5    = popt5[0]
psol_pow5 = popt5[1]
rmr_pow5  = popt5[2]
const6    = popt6[0]
psol_pow6 = popt6[1]
lne_pow6  = popt6[2]
rmr_pow6  = popt6[3]
const7    = popt7[0]
psol_pow7 = popt7[1]
lcon_pow7 = popt7[2]
const8 = popt8[0]; psol_pow8 = popt8[1]; lcon_pow8 = popt8[2]; rmr_pow8 = popt8[3]
const9 = popt9[0]; lcon_pow9 = popt9[1]
const10 = popt10[0]; psol_pow10 = popt10[1]; q95_pow10 = popt10[2]

# Variables for plotting.
otf_tot_fit1 = const1 * rmr ** rmr_pow1
otf_tot_fit2 = const2 * psol ** psol_pow2 * tri ** tri_pow2
otf_tot_fit3 = const3 * psol ** psol_pow3 * tri ** tri_pow3 * rmr ** rmr_pow3
otf_tot_fit4 = const4 * psol ** psol_pow4
otf_tot_fit5 = const5 * psol ** psol_pow5 * rmr ** rmr_pow5
otf_tot_fit6 = const6 * psol ** psol_pow6 * lne ** lne_pow6 * rmr ** rmr_pow6
otf_tot_fit7 = const7 * psol ** psol_pow7 * lcon ** lcon_pow7
otf_tot_fit8 = const8 * psol ** psol_pow8 * lcon ** lcon_pow8 * rmr ** rmr_pow8
otf_tot_fit9 = const9 * lcon ** lcon_pow9
otf_tot_fit10 = const10 * psol ** psol_pow10 * q95 ** q95_pow10
line_x = [0, 1]
line_y = [0, 1]

# Plot the scaling law results.
ms = 13
fig, ax = plt.subplots()
ax.plot(line_x, line_y, 'k--', lw=2)

# Choose correct points to plot.
if plot_scaling == 1:
    otf_tot_fit = otf_tot_fit1
    ax.set_xlabel(r"{:.2f}*$\mathdefault{{RMRSOMP^{{{:.1f}}}}}$".format(const1, rmr_pow1), fontsize=16)
elif plot_scaling == 2:
    otf_tot_fit = otf_tot_fit2
    ax.set_xlabel(r"{:.1f}*$\mathdefault{{PSOL^{{{:.1f}}}*TRITOP^{{{:.1f}}}}}$".format(const2, psol_pow2, tri_pow2), fontsize=16)
elif plot_scaling == 3:
    otf_tot_fit = otf_tot_fit3
    ax.set_xlabel(r"{:.1f}*$\mathdefault{{PSOL^{{{:.1f}}}*TRITOP^{{{:.1f}}}*RMRSOMP^{{{:.1f}}}}}$".format(const3, psol_pow3, tri_pow3, rmr_pow3), fontsize=16)
elif plot_scaling == 4:
    otf_tot_fit = otf_tot_fit4
    ax.set_xlabel(r"{:.2f}*$\mathdefault{{PSOL^{{{:.1f}}}}}$".format(const4, psol_pow4), fontsize=16)
elif plot_scaling == 5:
    otf_tot_fit = otf_tot_fit5
    ax.set_xlabel(r"{:.1f}*$\mathdefault{{PSOL^{{{:.1f}}}*RMRSOMP^{{{:.1f}}}}}$".format(const5, psol_pow5, rmr_pow5), fontsize=16)
elif plot_scaling == 6:
    otf_tot_fit = otf_tot_fit6
    ax.set_xlabel(r"{:.1f}*$\mathdefault{{PSOL^{{{:.1f}}}*\lambda^{{{:.1f}}}*RMRSOMP^{{{:.1f}}}}}$".format(const6, psol_pow6, lne_pow6, rmr_pow6), fontsize=16)
elif plot_scaling == 7:
    otf_tot_fit = otf_tot_fit7
    ax.set_xlabel(r"{:.2e} * $\mathdefault{{P_{{SOL}}^{{{:.1f}}} * L_{{CONN}}^{{{:.1f}}}}}$".format(const7, psol_pow7, lcon_pow7), fontsize=16)
elif plot_scaling == 8:
    otf_tot_fit = otf_tot_fit8
    ax.set_xlabel(r"{:.2e}*$\mathdefault{{PSOL^{{{:.1f}}}*LCONN^{{{:.1f}}}*RMRSOMP^{{{:.1f}}}}}$".format(const8, psol_pow8, lcon_pow8, rmr_pow8), fontsize=16)
elif plot_scaling == 9:
    otf_tot_fit = otf_tot_fit9
    ax.set_xlabel(r"{:.2e}*$\mathdefault{{LCONN^{{{:.1f}}}}}$".format(const9, lcon_pow9), fontsize=16)
elif plot_scaling == 10:
    otf_tot_fit = otf_tot_fit10
    ax.set_xlabel(r"{:.2e} * $\mathdefault{{P_{{SOL}}^{{{:.1f}}} * Q95^{{{:.1f}}}}}$".format(const10, psol_pow10, q95_pow10), fontsize=16)

ax.plot(otf_tot_fit[floor_mask], otf_tot[floor_mask], marker='*', lw=0, color="tab:purple", ms=ms, mec="k", label="Floor")
ax.plot(otf_tot_fit[shelf_mask], otf_tot[shelf_mask], marker='*', lw=0, color="tab:red", ms=ms, mec="k", label="Shelf")

# Plot where we would expect the not yet RBS'd probes to be.
if plot_none:
    if plot_scaling == 1:
        pass
    elif plot_scaling == 2:
        pass
    elif plot_scaling == 3:
        x = scaling3((psol_none, tri_none, rmr_none), *popt3)
        ax.scatter(x[floor_none], x[floor_none], s=200, marker="o", lw=3, facecolors='none', edgecolors="tab:purple")
        ax.scatter(x[shelf_none], x[shelf_none], s=200, marker="o", lw=3, facecolors='none', edgecolors="tab:red")
    elif plot_scaling == 4:
        pass

if labels:
    for i, txt in enumerate(probes):
        ax.annotate(txt, (otf_tot_fit[i], otf_tot[i]))
        #if plot_scaling == 1:
        #    ax.annotate(txt, (otf_tot_fit1[i], otf_tot[i]))
        #elif plot_scaling == 2:
        #    ax.annotate(txt, (otf_tot_fit2[i], otf_tot[i]))
        #elif plot_scaling == 3:
        #    ax.annotate(txt, (otf_tot_fit3[i], otf_tot[i]))
        #elif plot_scaling == 4:
        #    ax.annotate(txt, (otf_tot_fit4[i], otf_tot[i]))
    for i, txt in enumerate(probes_none):
        if plot_scaling == 1:
            pass
        elif plot_scaling == 2:
            pass
        elif plot_scaling == 3:
            ax.annotate(txt, (x[i], x[i]))
        elif plot_scaling == 4:
            pass

# Calculate R^2. Would just be an R^2 to the [0, 1] line.
y_true = otf_tot_fit
y_pred = otf_tot
#y_pred = otf_tot_fit
#sy_true = otf_tot
r2 = r2_score(y_true, y_pred)
print("R^2 = {:.3f}".format(r2))

ax.set_xlim([0, 0.4])
ax.set_ylim([0, 0.4])
ax.set_xticks(np.arange(0, 0.5, 0.1))
ax.set_yticks(np.arange(0, 0.5, 0.1))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(fontsize=14)
ax.set_ylabel("OTF Total W / Shot", fontsize=16)
ax.set_title(r"$\mathdefault{Bx \nabla B\downarrow}$ H-mode probes", fontsize=16)
fig.tight_layout()
fig.show()
