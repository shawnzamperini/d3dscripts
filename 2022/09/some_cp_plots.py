import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import numpy as np


root = "/Users/zamperini/My Drive/Research/Data/cp_data/"

# Scan through each file, load data accordingly.
cps = {}
for i in range(1, 11):
    if i < 10:
        numstr = "0{}".format(i)
    else:
        numstr = str(i)
    rpath = "{}MCPR{}W.csv".format(root, numstr)
    lpath = "{}MCPL{}W.csv".format(root, numstr)
    rdf = pd.read_csv(rpath)
    ldf = pd.read_csv(lpath)
    cps["MCPR{}W".format(i)] = rdf
    cps["MCPL{}W".format(i)] = ldf

extras = {}
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(9,4))
for i in range(1, 11):

    # Just grabbing the first shot for each R-Rsep values.
    if i == 1:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190422"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190422"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 2.2
        dens = 3.14e19

    elif i == 2:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190424"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190424"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 2.2
        dens = 3.47e19

    elif i == 3:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190426"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190426"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 2.2
        dens = 3.73e19

    elif i == 4:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190428"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190428"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 7.4
        dens = 6.50e19

    elif i == 5:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190450"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190450"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 6.8
        dens = 5.70e19

    elif i == 6:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190454"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190454"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 1.25

    elif i == 7:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190456"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190456"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 1.22
        dens = 3.45e19

    elif i == 8:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190459"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190459"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 1.22
        dens = 3.70e19

    elif i == 9:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190481"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190481"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 1.3
        dens = 3.7e19

    elif i == 10:
        xr = cps["MCPR{}W".format(i)]["R-Rsep OMP (cm) Shot:190491"]
        yr = cps["MCPR{}W".format(i)]["W areal density (1e15cm-2)"]
        xl = cps["MCPL{}W".format(i)]["R-Rsep OMP (cm) Shot:190491"]
        yl = cps["MCPL{}W".format(i)]["W areal density (1e15cm-2)"]
        pinj = 2.2
        dens = 3.7e19

    yr = medfilt(yr, 51)
    yl = medfilt(yl, 51)

    # Assign based off ITF/OTF.
    if i in [1,2,3,4,10]:
        xitf = xl
        xotf = xr
        yitf = yl
        yotf = yr
    else:
        xitf = xr
        xotf = xl
        yitf = yr
        yotf = yl

    # Assign some additional data.
    extras[i] = {}
    extras[i]["maxwitf"] = yitf.max()
    extras[i]["maxwotf"] = yotf.max()
    extras[i]["pinj"] = pinj
    extras[i]["dens"] = dens

    # The values at a point away from re-erosion.
    idx1 = np.argmin(np.abs(xitf-12))
    idx2 = np.argmin(np.abs(xotf-12))
    extras[i]["winditf"] = yitf[idx1]
    extras[i]["windotf"] = yotf[idx2]

    label = i
    ax1.plot(xitf, yitf, label=label)
    ax2.plot(xotf, yotf, label=label)

ax2.legend()
fig.supxlabel("R-Rsep OMP (cm)", fontsize=14)
ax1.set_ylabel("W Areal Density (W/cm2)", fontsize=14)
ax1.set_title("ITF", fontsize=14)
ax2.set_title("OTF", fontsize=14)
fig.tight_layout()
fig.show()

# Additional plots.
maxitfotfs = []; winditfotfs = []
pinjs = []; denss = []
for i in range(1, 11):
    maxitfotfs.append(extras[i]["maxwitf"] / extras[i]["maxwotf"])
    winditfotfs.append(extras[i]["winditf"] / extras[i]["windotf"])
    pinjs.append(extras[i]["pinj"])
    denss.append(extras[i]["dens"])

fig, ax1 = plt.subplots(figsize=(5,4))
#ax1.scatter(pinjs, maxitfotfs)

# Excluding MCP10 because it's very odd.
for i in range(1, 10):

    # Forward (unfavorable)
    if i in [1,2,3,4,10]:
        if i == 1:
            ax1.scatter(denss[i-1], winditfotfs[i-1], marker="o", color="tab:purple", label="Unfavorable")
        else:
            ax1.scatter(denss[i-1], winditfotfs[i-1], marker="o", color="tab:purple")
        #ax1.scatter(pinjs[i-1], winditfotfs[i-1], marker="o")
    else:
        if i == 5:
            ax1.scatter(denss[i-1], winditfotfs[i-1], marker="o", color="tab:red", label="Favorable")
        else:
            ax1.scatter(denss[i-1], winditfotfs[i-1], marker="o", color="tab:red")
        #ax1.scatter(pinjs[i-1], winditfotfs[i-1], marker=".")
#ax1.set_xlabel("PINJ (MW)", fontsize=14)
ax1.legend(fontsize=14)
ax1.set_xlabel("Line-averaged density (m-3)", fontsize=14)
ax1.set_ylabel("ITF/OTF (Windowed region)", fontsize=14)
fig.tight_layout()
fig.show()
