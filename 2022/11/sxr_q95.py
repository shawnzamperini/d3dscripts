# Is there any trend in the Fe23 signal from SPRED with q95 (indicative of
# tile misalignments)?
from gadata import gadata
import numpy as np
import matplotlib.pyplot as plt
import MDSplus
from scipy.interpolate import interp1d
from scipy.signal import medfilt


shots = [190446, 190447, 190450, 190451, 190428, 190429, 190055, 190056, 190059,
    190061, 190996, 190999, 191001, 191002]
# 190502: Taken out because ELM surpression lead to high SXR levels.
# 190448: Density scan
# 191014: Rough density scan
# Forward BT: 190428 190429
ech = [190056, 190059, 190061, 190996, 190999, 191001, 191002]

# Time range once strike point is parked on progressive angle.
tmin = 2600
tmax = 3500

results = {}
conn = MDSplus.Connection("atlas.gat.com")
for shot in shots:

    # Load signals.
    print(shot)
    gaobj = gadata("SX195R1S2", shot, connection=conn)
    sxr_t = np.array(gaobj.xdata)
    sxr = np.array(gaobj.zdata)
    gaobj = gadata("Q95", shot, connection=conn)
    q95_t = gaobj.xdata
    q95 = gaobj.zdata

    # Need to smooth fe23.
    sxrs = medfilt(sxr, 21)

    # Need to interpolate onto a common time axis.
    f_q95 = interp1d(q95_t, q95)

    # Get values of each just for the time range.
    sxr_mask = np.logical_and(sxr_t>tmin, sxr_t<tmax)
    sxr_keep = sxr[sxr_mask]
    sxrs_keep = sxrs[sxr_mask]
    q95_keep = f_q95(sxr_t[sxr_mask])

    # Find max fe23 and the corresponding q95 value.
    max_idx = np.argmax(sxrs_keep)
    sxr_max = sxrs_keep[max_idx]
    q95_max = q95_keep[max_idx]

    # Grab the Langmuir probe data and pair each data point up with its
    # respective Te value at that time.
    gaobj = gadata(".PROBE_073.TEMP", shot, connection=conn, tree="LANGMUIR")
    lp_t = np.array(gaobj.xdata)
    lp_te = np.array(gaobj.zdata)
    lp_tes = medfilt(lp_te, 51)
    f_lp_te = interp1d(lp_t, lp_tes)
    lp_te_int = f_lp_te(sxr_t[sxr_mask])

    gaobj = gadata(".PROBE_073.DENS", shot, connection=conn, tree="LANGMUIR")
    lp_t = np.array(gaobj.xdata)
    lp_ne = np.array(gaobj.zdata)
    lp_nes = medfilt(lp_ne, 51)
    f_lp_ne = interp1d(lp_t, lp_nes)
    lp_ne_int = f_lp_ne(sxr_t[sxr_mask]) * 1e6  # cm-3 to m-3

    # Likewise with the line-averaged density.
    gaobj = gadata("DENSV2", shot, connection=conn)
    ne_t = np.array(gaobj.xdata)
    ne = np.array(gaobj.zdata)
    nes = medfilt(ne, 21)
    f_ne = interp1d(ne_t, nes)
    ne_int = f_ne(sxr_t[sxr_mask])

    # Change in stored energy.
    if shot in [190428, 190996, 190999, 191001, 191002]:
        mhd_int = np.zeros(len(sxr_t[sxr_mask]))
    else:
        gaobj = gadata("DWMHDF", shot, connection=conn)
        mhd_t = np.array(gaobj.xdata)
        mhd = np.array(gaobj.zdata)
        f_mhd = interp1d(mhd_t, mhd)
        mhd_int = f_mhd(sxr_t[sxr_mask])

    # Filterscope data.
    gaobj = gadata("TUBE38:PMT_VOLT", shot, connection=conn, tree="FSCOPE")
    fs_t = np.array(gaobj.xdata)
    fs = np.array(gaobj.zdata)
    f_fs = interp1d(fs_t, fs)
    fs_int = f_fs(sxr_t[sxr_mask])

    # Separatrix Te.
    gaobj = gadata("ZUPERTS", shot, connection=conn)
    tsz_t = np.array(gaobj.xdata)
    tsz = np.array(gaobj.zdata)
    tsz_avg = tsz[np.logical_and(tsz_t>tmin, tsz_t<tmax)].mean()

    # Load closest TS chord to average separatrix location. Drop zeros.
    if shot not in [190428]:
        gaobj = gadata("TSTE_{:}".format(tsz_avg*100), shot, connection=conn)
        tste_t = np.array(gaobj.xdata)
        tste = np.array(gaobj.zdata)
        tmpmask = np.nonzero(tste)
        tste = tste[tmpmask]
        tste_t = tste_t[tmpmask]
        f_tste = interp1d(tste_t, tste)
        tste_int = f_tste(sxr_t[sxr_mask])
    else:
        tste_int = np.zeros(len(sxr_t[sxr_mask]))

    # Particle confinement time.
    gaobj = gadata("TAUMHD", shot, connection=conn)
    tau_t = np.array(gaobj.xdata)
    tau = np.array(gaobj.zdata)
    taus = medfilt(tau, 21)
    f_tau = interp1d(tau_t, taus)
    tau_int = f_tau(sxr_t[sxr_mask])


    results[shot] = {"q95":q95_keep, "sxr":sxr_keep, "sxrs":sxrs_keep,
        "sxr_max":sxr_max, "q95_max":q95_max, "t":sxr_t[sxr_mask],
        "max_idx":max_idx, "lp_t":lp_t, "lp_te":lp_te, "lp_tes":lp_tes,
        "lp_te_int":lp_te_int, "ne_t":ne_t, "ne":ne, "nes":nes, "ne_int":ne_int,
        "lp_ne_int":lp_ne_int, "mhd_int":mhd_int, "fs_int":fs_int, "tste_int":tste_int,
        "tau_int":tau_int}

fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, figsize=(12,7))
colors = ["C{:}".format(i) for i in range(0, len(shots))]
count = 0
for shot, data in results.items():

    # Deuterium flux to target so we can normalize.
    mi = 931.49e6
    dflux = data["lp_ne_int"] * np.sqrt(2*data["lp_te_int"]/mi) * 3e8

    ax1.plot(data["t"], data["sxr"], alpha=0.3, color=colors[count], zorder=5)
    ax1.plot(data["t"], data["sxrs"], color=colors[count], zorder=10)
    #ax1.scatter(data["t"][data["max_idx"]], data["sxrs"][data["max_idx"]], color=colors[count], marker="*", edgecolor="k", zorder=15)
    ax2.plot(data["t"], data["tste_int"], color=colors[count], label=shot)
    #ax3.scatter(data["q95_max"], data["sxr_max"], label=shot, color=colors[count])
    ax3.scatter(data["q95"], data["sxr"]/dflux, color=colors[count], label=shot, s=1, alpha=0.3)

    ax4.plot(data["lp_t"], data["lp_te"], alpha=0.3, color=colors[count], zorder=5)
    ax4.plot(data["lp_t"], data["lp_tes"], color=colors[count], zorder=10)

    ax5.scatter(data["lp_te_int"], data["sxr"], color=colors[count], label=shot, s=1, alpha=0.3)
    #ax5.scatter(data["lp_te_int"].mean(), data["sxr"].mean(), color=colors[count], label=shot, s=15)
    #ax5.scatter(data["lp_te_int"], data["sxr"]/data["tau_int"], color=colors[count], label=shot, s=1, alpha=0.3)

    #if shot not in [190428]:

    #    if shot in ech:
    #        ax5.scatter(data["lp_ne_int"]*data["lp_te_int"], data["sxr"]/(data["tste_int"]/data["lp_te_int"]).mean(), color=colors[count], label=shot, s=1, alpha=0.3)
    #    else:
    #        ax5.scatter(data["lp_ne_int"]*data["lp_te_int"], data["sxr"]/(data["tste_int"]/data["lp_te_int"]).mean(), color=colors[count], label=shot, s=1, alpha=0.3)

    #if shot not in [190428, 190996, 190999, 191001, 191002]:
    #    ax5.scatter(data["mhd_int"].mean(), data["sxr"].mean()/dflux.mean(), color=colors[count], label=shot, s=15)

    ax6.scatter(data["q95"], data["lp_te_int"], color=colors[count], s=1, alpha=0.3)

    ax7.plot(data["ne_t"], data["ne"], alpha=0.3, color=colors[count], zorder=5)
    ax7.plot(data["ne_t"], data["nes"], color=colors[count], zorder=10)
    ax8.scatter(data["fs_int"].mean(), data["sxr"].mean(), color=colors[count], label=shot, s=15)

    count += 1

#ax5.set_ylim([0, 1])
ax2.legend()
ax3.set_xlabel("Te Sep (eV)")
ax3.set_ylabel("SXR Signal")
ax5.set_xlabel("Te Targ")
ax5.set_ylabel("SXR Signal")
ax6.set_xlabel("q95")
ax6.set_ylabel("Te Targ")
ax8.set_xlabel("VW2")
ax8.set_ylabel("SXR")
fig.tight_layout()
fig.show()


# Now let's just run a simple linear regression on the data to see if it aligns.
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split


# Create a number of train/test datasets where for each shot is a test dataset,
# and the training set is the rest of the shots. We train the same number of
# models and combine the results.
ypreds = np.array([])
ytests = np.array([])
for testshot, testdata in results.items():
    print(testshot)

    # Assemble training set out of all the other shots.
    sxr = np.array([])
    ne = np.array([])
    lpte = np.array([])
    lpne = np.array([])
    tste = np.array([])
    q95 = np.array([])
    tau = np.array([])
    for trainshot, traindata in results.items():
        if trainshot == testshot:
            continue
        sxr = np.append(sxr, traindata["sxrs"])
        ne = np.append(ne, traindata["ne_int"])
        lpte = np.append(lpte, traindata["lp_te_int"])
        lpne = np.append(lpne, traindata["lp_ne_int"])
        tste = np.append(tste, traindata["tste_int"])
        q95 = np.append(q95, traindata["q95"])
        tau = np.append(tau, traindata["tau_int"])
    leak_proxy = tste / lpte
    y_train = sxr
    #X_train = X = np.vstack([ne, lpte, q95]).T
    #X_train = X = np.vstack([lpte, lpne, leak_proxy]).T
    X_train = X = np.vstack([lpte, tau]).T
    y_test = testdata["sxr"]
    #X_test = np.vstack([testdata["ne_int"], testdata["lp_te_int"], testdata["q95"]]).T
    #X_test = np.vstack([testdata["lp_te_int"], testdata["lp_ne_int"], testdata["tste_int"]/testdata["lp_te_int"]]).T
    X_test = np.vstack([testdata["lp_te_int"], testdata["tau_int"]]).T

    # Create model and predict.
    scaler = StandardScaler().fit(X_train)
    X_train_scaled = scaler.transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    #model = LinearRegression().fit(X_train_scaled, y_train)
    print("  training model...")
    model = RandomForestRegressor(n_jobs=-1).fit(X_train_scaled, y_train)
    ypred = model.predict(X_test_scaled)

    ypreds = np.append(ypreds, ypred)
    ytests = np.append(ytests, y_test)

fig, ax = plt.subplots()
ax.plot([0, 2], [0, 2], color="k", linestyle="--")
ax.scatter(ypreds, ytests, s=5, alpha=0.25, color="k")
ax.set_xlabel("model")
ax.set_ylabel("actual")
ax.set_xlim([0, 1.7])
ax.set_ylim([0, 1.7])
ax.grid()
fig.tight_layout()
fig.show()


"""
# Consolidate data.
sxr = np.array([])
ne = np.array([])
lpte = np.array([])
q95 = np.array([])
for shot, data in results.items():
    sxr = np.append(sxr, data["sxrs"])
    ne = np.append(ne, data["ne_int"])
    lpte = np.append(lpte, data["lp_te_int"])
    q95 = np.append(q95, data["q95"])
y = sxr
#X = np.vstack([q95]).T
#X = np.vstack([lpte, q95]).T
X = np.vstack([ne, lpte, q95]).T
X_train, X_test, y_train, y_test = train_test_split(X, y)

# Create model.
scaler = StandardScaler().fit(X)
X_train_scaled = scaler.transform(X_train)
X_test_scaled = scaler.transform(X_test)
#model = LinearRegression().fit(X_train_scaled, y_train)
model = RandomForestRegressor().fit(X_train_scaled, y_train)

# See how well it does.
ypred = model.predict(X_test_scaled)

fig, ax = plt.subplots()
ax.plot([0, 2], [0, 2], color="k", linestyle="--")
ax.scatter(ypred, y_test, s=5, color="k")
ax.set_xlabel("model")
ax.set_ylabel("actual")
ax.set_xlim([0, 1.7])
ax.set_ylim([0, 1.7])
ax.grid()
fig.tight_layout()
fig.show()
"""
