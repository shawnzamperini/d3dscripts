# Script to make some nice plots with errorbars of the far-SOL RCP data,
# searching for trends in the falloff width, and maybe even Iwall.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression


# Constants.
root = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/"
shots = [190440, 190442, 190484, 190485, 190486]
fit_region = [228.9, 232.0]

def exp_fit(x, a, b, c):
    return a * np.exp(-x*b) + c

data = []
for shot in shots:
    for i in [1, 2]:

        # Load, extract data.
        df = pd.read_csv(root+"MP{}_{}.tab".format(shot, i), delimiter="\t")
        r = df["R(cm)"]
        ne = df["Ne(E18 m-3)"]
        te = df["Te(eV)"]
        m = df["Machn"]
        shot_label = "{} #{}".format(shot, i)
        print(shot_label)

        # Assign shot parameters. Some shots may require special hands on fixes.
        offset_flag = False
        if shot == 190440 and i == 1:
            offset = np.abs(ne[8])
            offset_flag = True
            nebar = 1.88e19
            vblob = 473
            mach = 0.0611
            zvsout = 1.19242
            jsat = 0.232
            fblob = 1.74e4
            dgap = 0.015
        elif shot == 190440 and i == 2:
            nebar = 2.74e19
            vblob = 286
            mach = 0.0478
            zvsout = 1.21762
            jsat = 0.204
            fblob = 1.10e5
            dgap = -0.015
        elif shot == 190442 and i == 1:
            nebar = 1.64e19
            vblob = 194
            mach = -0.3373
            zvsout = 1.22096
            jsat = 0.162
            fblob = 3.81e4
            dgap = 0.038
        elif shot == 190442 and i == 2:
            nebar = 3.66e19
            vblob = 373
            mach = -0.3026
            zvsout = 1.22114
            jsat = 0.716
            fblob = 5.5e4
            dgap = -0.038
        elif shot == 190484 and i == 1:
            offset = np.abs(ne[9])
            offset_flag = True
            nebar = 2.97e19
            vblob = 590
            mach = 0.2220
            zvsout = 1.21986
            jsat = 0.265
            fblob = 9.9e4
            dgap = 0.004
        elif shot == 190484 and i == 2:
            offset = np.abs(ne[8])
            offset_flag = True
            nebar = 2.97e19
            vblob = 327
            mach = 0.1351
            zvsout = 1.22301
            jsat = 0.276
            fblob = 4.71e4
            dgap = -0.004
        elif shot == 190485 and i == 1:
            nebar = 3.47e19
            vblob = 448
            mach = -0.3940
            zvsout = 1.22256
            jsat = 0.413
            fblob = 3.3e4
            dgap = 0.008
        elif shot == 190485 and i == 2:
            nebar = 3.47e19
            vblob = 469
            mach = -0.4274
            zvsout = 1.22413
            jsat = 0.526
            fblob = 9.9e4
            dgap = -0.008
        elif shot == 190486 and i == 1:
            nebar = 4.03e19
            vblob = 324
            mach = 0.0223
            zvsout = 1.22150
            jsat = 0.501
            fblob = 2.61e4
            dgap = -0.009
        elif shot == 190486 and i == 2:
            nebar = 4.03e19
            vblob = 401
            mach = 0.0253
            zvsout = 1.22371
            jsat = 0.572
            fblob = 1.1e5
            dgap = -0.009

        if offset_flag:
            print("{} #{}".format(shot, i))
            print(" Non-zero offset applied: {:.2}".format(offset))
            ne += offset

        # Exponential fit in the fit region. To get some form of statistics,
        # we will fit multiple regions defined by  +/- a set amount (permutate
        # that amount to get multiple windows).
        r_width = 0.0
        r_help = 229
        mask1 = np.logical_and(r>=fit_region[0], r<=fit_region[1])
        mask2 = np.logical_and(r>=fit_region[0]-r_width, r<=fit_region[1])
        mask3 = np.logical_and(r>=fit_region[0], r<=fit_region[1]-r_width)
        mask4 = np.logical_and(r>=fit_region[0]+r_width, r<=fit_region[1])
        mask5 = np.logical_and(r>=fit_region[0]+r_width, r<=fit_region[1]+r_width)
        mask6 = np.logical_and(r>=fit_region[0]-r_width, r<=fit_region[1]+r_width)
        mask7 = np.logical_and(r>=fit_region[0]+r_width, r<=fit_region[1]-r_width)
        mask8 = np.logical_and(r>=fit_region[0], r<=fit_region[1]+r_width)
        mask9 = np.logical_and(r>=fit_region[0]-r_width, r<=fit_region[1]-r_width)
        masks = [mask1, mask2, mask3, mask4, mask5, mask6, mask7, mask8, mask9]
        lamb_fars = []; popts = []; ints_far = []; models = []
        count = 1
        for mask in masks:
            #try:

                # Fitting an exponential.
                #popt, pcov = curve_fit(exp_fit, r[mask] - r_help, ne[mask], p0=(1, -20, 0))
                #rfit = np.linspace(fit_region[0]-r_width, fit_region[1]+r_width, 25)
                #nefit = exp_fit(rfit-r_help, *popt)
                #lamb_fars.append(-1 / popt[1])
                #popts.append(popt)

            # Or fitting a line to the ln(ne).
            model = LinearRegression().fit(r[mask].values.reshape((-1, 1)), np.log(ne[mask]).values)
            lamb_fars.append(-1/model.coef_[0])
            ints_far.append(model.intercept_)
            models.append(model)

            #except:
            #    print("{}: Unable to fit mask{}".format(shot_label, count))
            count += 1
        int_far_mean = np.mean(ints_far)
        lamb_far_mean = np.mean(lamb_fars)
        lamb_far_std = np.std(lamb_fars)

        data.append({"shot_label":shot_label, "r":r, "ne":ne, "te":te,
            "nebar":nebar, "m":m, "lamb_fars":lamb_fars,
            "lamb_far_mean":lamb_far_mean, "lamb_far_std":lamb_far_std,
            "vblob":vblob, "mach":mach, "zvsout":zvsout, "jsat":jsat,
            "int_far_mean":int_far_mean, "models":models, "fblob":fblob,
            "dgap":dgap})

# Assign colors based off nebar.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
nebars = []; nebar_min = 1e25; nebar_max = 0
for i in range(0, len(data)):
    nebar = data[i]["nebar"]
    if nebar < nebar_min:
        nebar_min = nebar
    if nebar > nebar_max:
        nebar_max = nebar

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))
ax1.axhline(0, color="k")
for i in range(0, len(data)):
    r = data[i]["r"]
    ne = data[i]["ne"]
    nebar = data[i]["nebar"]
    shot_label = data[i]["shot_label"]
    vblob = data[i]["vblob"]
    lamb_mean = data[i]["lamb_far_mean"]
    lamb_std = data[i]["lamb_far_std"]
    mach = data[i]["mach"]
    zvsout = data[i]["zvsout"]
    jsat = data[i]["jsat"]
    fblob = data[i]["fblob"]
    dgap = data[i]["dgap"]
    c = cmap((nebar - nebar_min)/(nebar_max - nebar_min))
    #ax1.scatter(r, ne*1e18, color=c, label=shot_label, edgecolors="k", marker="^")
    ax1.scatter(r, ne*1e18, label=shot_label, edgecolors="k", marker="^")
    ax2.errorbar(fblob, lamb_mean, yerr=lamb_std, marker="^", ms=10, lw=0, mec="k", elinewidth=2)
    if shot_label[-1] == "2":
        lamb_change = lamb_mean - data[i-1]["lamb_far_mean"]
        ax3.errorbar(dgap, lamb_change, yerr=lamb_std, marker="^", ms=10, lw=0, mec="k", elinewidth=2, label=shot_label[:-2])

ax1.set_xlabel("R (cm)")
ax1.set_ylabel("ne (m-3)")
ax1.legend()
ax1.set_ylim([1e17, 1.5e19])
ax1.grid(which="both", alpha=0.3)
ax1.set_yscale("log")
ax2.grid()
ax2.set_ylabel("lambda far-SOL (cm)")
ax3.legend()

fig.tight_layout()
fig.show()


# Just a temporary plot for 190485 (i = 6 and 7).
fig, ax1 = plt.subplots(figsize=(5,4))
idx1 = 6
idx2 = 7
x1 = data[idx1]["r"]
y1 = data[idx1]["ne"]
x2 = data[idx2]["r"]
y2 = data[idx2]["ne"]
l1 = data[idx1]["lamb_far_mean"]
l2 = data[idx2]["lamb_far_mean"]
i1 = data[idx1]["int_far_mean"]
i2 = data[idx2]["int_far_mean"]
model1 = data[idx1]["models"][0]
model2 = data[idx2]["models"][0]
rfit = np.linspace(228.9, 232, 10)
nefit1 = np.exp(model1.predict(rfit.reshape(-1, 1)))
nefit2 = np.exp(model2.predict(rfit.reshape(-1, 1)))
ax1.axvline(228.9, color="k", linestyle="--")
ax1.scatter(x1, y1*1e18, label="3000 ms", edgecolors="k", marker="^", color="tab:red")
ax1.scatter(x2, y2*1e18, label="4000 ms", edgecolors="k", marker="^", color="tab:purple")
ax1.plot(rfit, nefit1*1e18, color="tab:red")
ax1.plot(rfit, nefit2*1e18, color="tab:purple")
ax1.set_xlabel("R (cm)")
ax1.set_ylabel("ne (m-3)")
ax1.legend()
ax1.set_ylim([3e17, 1.5e19])
ax1.grid(which="both", alpha=0.3)
ax1.set_yscale("log")
fig.tight_layout()
fig.show()
