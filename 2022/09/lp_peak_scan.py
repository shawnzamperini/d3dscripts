# The gist here is to get the disatnce from the peak LP ne value to the strike
# point as a function of neSOL.
import numpy as np
import matplotlib.pyplot as plt
import get_lp
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy import interpolate
from scipy.signal import medfilt


# Candidates
# 167195 4000-5000

shot = 174175
tstart = 2000
tend = 3600
region = "shelf"
fit_range = [0.0, 0.2]
lpdict = get_lp.plot_lps(shot, tstart, tend, bins=100, tunnel=False, showplot=False)


# Plot only SAS probes.
if region == "sas":
    labels = ["A-4", "A-5", "A-6", "A-7", "A-8", "A-9", "A10", "A11", "A12",
        "A13", "A14", "A20", "A21", "A22"]
elif region == "shelf":
    labels = ["S-1", "S-2", "S-3", "S-4", "S-5", "S-6", "S-7", "S-8", "S-9",
        "S10", "S11"]

colors = {}
for i in range(0, len(labels)):
    colors[labels[i]] = "C{}".format(i)

data = {"rmrs":[], "ne":[], "te":[], "color":[]}
for i in range(0, len(lpdict["labels"])):

    label = lpdict["labels"][i].strip()
    if label in labels:
        data["rmrs"].append(lpdict["rminrsep"][i])
        data["ne"].append(lpdict["ne (cm-3)"][i] * 1e6)
        data["te"].append(lpdict["Te (eV)"][i])
        data["color"].append(colors[label])

sort_idx = np.argsort(data["rmrs"])
data["rmrs"] = np.array(data["rmrs"])[sort_idx]
data["ne"] = np.array(data["ne"])[sort_idx]
data["te"] = np.array(data["te"])[sort_idx]
data["color"] = np.array(data["color"])[sort_idx]

# Fit a conv guass to the region to find the location of the peak value.
def gauss_conv_exp_fit(s, width, lambda_n, n0, n_bg, s0):
    fx = 5
    return n0 / 2.0 * np.exp((width/(2*lambda_n*fx))**2 - (s-s0)/(lambda_n *
      fx)) * erfc(width/(2*lambda_n*fx) - (s-s0)/width)

fit_idx = np.logical_and(data["rmrs"]>=fit_range[0], data["rmrs"]<=fit_range[1])
#guess = (0.05, 0.02, 1.0, 0.0, 0.0)
#popt, pcov = curve_fit(gauss_conv_exp_fit, data["rmrs"][fit_idx],
#  data["ne"][fit_idx]/1e18, p0=guess, maxfev=5000)
#fit_rmrs = data["rmrs"][fit_idx]
#fit_ne = gauss_conv_exp_fit(fit_rmrs, *popt)

# Median filter to the data.
nefilt = medfilt(data["ne"], 35)
tefilt = medfilt(data["te"], 35)
data["nefilt"] = nefilt
data["tefilt"] = tefilt

# Peak ne occurs at...
peak_idx = np.argmax(nefilt)
rpeak = data["rmrs"][peak_idx]
nepeak = nefilt[peak_idx]
print("Peak value at {:.2f} cm from the strike point".format(rpeak*100))

# Spline interpolation, but I don't understand it too well.
knot_numbers = 5
x_new = np.linspace(0, 1, knot_numbers+2)[1:-1]
q_knots = np.quantile(data["rmrs"], x_new)
t, c, k = interpolate.splrep(data["rmrs"], nefilt/1e18, t=q_knots, s=1)
nefit = interpolate.BSpline(t, c, k)(data["rmrs"]) * 1e18

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

ax1.scatter(data["rmrs"], data["ne"], c=data["color"], zorder=5)
ax1.plot(data["rmrs"], nefilt, color="k", zorder=10)
#ax1.scatter(rpeak, nepeak, marker="*", s=200, color="r", edgecolors="k", zorder=15)
#ax1.plot(data["rmrs"], nefit, color="r")
ax1.set_xlabel("R-Rsep (m)", fontsize=14)
ax1.set_ylabel("ne (m-3)", fontsize=14)

ax2.scatter(data["rmrs"], data["te"], c=data["color"])
ax2.plot(data["rmrs"], tefilt, color="k")
ax2.set_xlabel("R-Rsep (m)", fontsize=14)
ax2.set_ylabel("Te (eV)", fontsize=14)

fig.tight_layout()
fig.show()
