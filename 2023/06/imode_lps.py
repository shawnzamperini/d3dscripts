import get_lp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Values of when I-mode occurred according to Amanda.
imode = {189384: [2700, 4000],
         189376: [[3400, 3750], [2800, 3100]],
         189377: [2800, 3200],
         189379: [2800, 3200],
         189380: [[3400, 3700], [2600, 3000]],
         189381: [2700, 3000],
         189382: [3100, 4300],
         189383: [3100, 3300]}


# Load data for each configuration.
path = "/Users/zamperini/My Drive/Research/Documents/2023/06/D3D I-mode shotlist 2022.xlsx"
df = pd.read_excel(path, skiprows=1)
lsn = df.iloc[0:14]
usn = df.iloc[15:18]

# Load each LP data.
lsn_dict = {}
lsn_lmode_dict = {}
usn_dict = {}
for i in range(0, len(lsn)):
    s = lsn.iloc[i]
    print(s["Shot"])
    lps = get_lp.plot_lps(s["Shot"], s["t1"]*1000, s["t1.1"]*1000, xtype="psin", showplot=False, tunnel=False)
    if s["Shot"] in lsn_dict.keys():
        lsn_dict[str(s["Shot"]) + " (1)"] = lsn_dict.pop(s["Shot"])
        lsn_dict[str(s["Shot"]) + " (2)"] = lps
    else:
        lsn_dict[s["Shot"]] = lps

    # Likewise but for the L-mode data beforehand (only for the shots where we have the associated heat flux data from S-4).
    if ~np.isnan(s["avg S-4 heat flux (W/cm2)"]):
        print("Loading L-mode data range...")
        lps = get_lp.plot_lps(s["Shot"], s["t2"] * 1000, s["t2.1"] * 1000, xtype="psin", showplot=False, tunnel=False)
        if s["Shot"] in lsn_lmode_dict.keys():
            lsn_lmode_dict[str(s["Shot"]) + " (1)"] = lsn_lmode_dict.pop(s["Shot"])
            lsn_lmode_dict[str(s["Shot"]) + " (2)"] = lps
        else:
            lsn_lmode_dict[s["Shot"]] = lps


for i in range(0, len(usn)):
    s = usn.iloc[i]
    lps = get_lp.plot_lps(s["Shot"], s["t1"]*1000, s["t1.1"]*1000, xtype="psin", showplot=False, tunnel=False)
    if s["Shot"] in usn_dict.keys():
        usn_dict[str(s["Shot"]) + " (1)"] = usn_dict.pop(s["Shot"])
        usn_dict[str(s["Shot"]) + " (2)"] = lps
    else:
        usn_dict[s["Shot"]] = lps


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 6))

# For the LSN let's plot all the shelf probes.
avg_heats = {}
avg_lmode_heats = {}
for shot, data in lsn_dict.items():
    psin = []
    te = []
    ne = []
    for i in range(0, len(data["psin"])):
        if data["labels"][i].strip()[0] == "S":
            psin.append(data["psin"][i])
            te.append(data["Te (eV)"][i])
            ne.append(data["ne (cm-3)"][i] * 1e6)

    ax1.scatter(psin, te, label=shot)
    ax3.scatter(psin, ne, label=shot)

    # Additional statistics.
    tmp_hf = []
    for i in range(0, len(data["psin"])):
        if data["labels"][i].strip() == "S-4":
            tmp_hf.append(data["heatflux (W/cm2)"][i])
    avg_heats[shot] = np.mean(tmp_hf)

# Just want the heat flux data from L-mode.
avg_lmode_heats = {}
for shot, data in lsn_lmode_dict.items():
    tmp_hf = []
    for i in range(0, len(data["psin"])):
        if data["labels"][i].strip() == "S-4":
            tmp_hf.append(data["heatflux (W/cm2)"][i])
    avg_lmode_heats[shot] = np.mean(tmp_hf)

# And for USN all the "E" and "P" probes.
for shot, data in usn_dict.items():
    psin = []
    te = []
    ne = []
    for i in range(0, len(data["psin"])):
        if data["labels"][i].strip()[0] in ["E", "P"]:
            psin.append(data["psin"][i])
            te.append(data["Te (eV)"][i])
            ne.append(data["ne (cm-3)"][i] * 1e6)

    ax2.scatter(psin, te, label=shot)
    ax4.scatter(psin, ne, label=shot)

ax1.set_ylabel("Te (eV)")
ax3.set_ylabel("ne (m-3)")
ax1.set_title("LSN")
ax2.set_title("USN")
ax1.legend()
ax2.legend()
fig.tight_layout()
fig.show()