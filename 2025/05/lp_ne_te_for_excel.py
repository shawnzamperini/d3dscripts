# Purpose of this script is to pull LP data for a given shot/time and then
# plot some various exponential curves on top to get a sense of what the range
# of falloff lengths (from the peak value) could be
import get_lp
import matplotlib.pyplot as plt
import numpy as np


# Inputs
#shot = 200083
#time = 2400
#time_window = 250
#peak_ne = 8e18
#peak_te = 50
shot = 201661
time = 3000
time_window = 250
peak_ne = 1.5e19
peak_te = 50
lambs = np.arange(0.04, 0.21, 0.02)

# Load LP data
lps = get_lp.plot_lps(shot, time - time_window / 2, time + time_window / 2, 
	tunnel=False, showplot=False)

"""
# Extract all the shelf probes (label starts with "S")
shelf_probes = []
for p, pdata in lp_dict.items():
	label = pdata["label"].strip()
	if label[0] == "S":
		shelf_probes.append(pdata)

# Extract data of dist_from_sep, ne and Te
drsep = []
ne = []
te = []
ps = []
for p in shelf_probes:

	# Limit to time window
	mask = np.logical_and(p["time"] > (time - time_window / 2), 
		p["time"] < time + time_window / 2)

	drsep = np.append(drsep, np.array(p["delrsepout"])[mask])
	ne = np.append(ne, np.array(p["dens"])[mask])
	te = np.append(te, np.array(p["temp"])[mask])
	ps = np.append(ps, np.full(mask.sum(), p["label"].strip()))
"""

# Select shelf probes
#mask = np.array([l.strip()[0]=="S" for l in lps["labels"]])
mask = np.array([l.strip()[0]=="F" for l in lps["labels"]])

# Get shelf data
rminrsep = np.array(lps["rminrsep"])[mask]
ne = np.array(lps["ne (cm-3)"])[mask] * 1e6  # cm-3 to m-3
te = np.array(lps["Te (eV)"])[mask]

# Generate some exponential data to plot over top for reference
xexp = np.linspace(-0.01, 0.20, 100)
yexp_ne = []
yexp_te = []
for l in lambs:
	yexp_ne.append(peak_ne * np.exp(-xexp / l))
	yexp_te.append(peak_te * np.exp(-xexp / l))



fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

ax1.set_xlabel("R-Rsep (m)")
ax2.set_xlabel("R-Rsep (m)")
ax1.set_ylabel("ne (m-3)")
ax2.set_ylabel("Te (eV)")

ax1.axvline(0.0, color="k")
ax1.scatter(rminrsep, ne)
for i in range(0, len(lambs)):
	ax1.plot(xexp, yexp_ne[i], label="{:.2f}".format(lambs[i]))

ax2.axvline(0.0, color="k")
ax2.scatter(rminrsep, te)
for i in range(0, len(lambs)):
	ax2.plot(xexp, yexp_te[i], label="{:.2f}".format(lambs[i]))

ax2.legend()
fig.tight_layout()
fig.show()
