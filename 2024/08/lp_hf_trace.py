import get_lp
import matplotlib.pyplot as plt
import numpy as np


shot = 200073
tmin = 3000
tmax = 3500

time = np.array([])
qpar = np.array([])
psin = np.array([])

# Not efficient, but run this X times to get different time ranges
tstep = 100
results = []
times = []
t = tmin
while (t < tmax):
    lpdict = get_lp.plot_lps(shot, t, t + tstep, xtype="psin", tunnel=False, showplot=False)
    results.append(lpdict)
    times.append(t + tstep / 2)
    t += tstep

    # Add into master arrays
    psin = np.append(psin, lpdict["psin"])
    qpar = np.append(qpar, lpdict["heatflux (W/cm2)"])
    time = np.append(time, np.full(len(lpdict["psin"]), t + tstep / 2))

fig, ax1 = plt.subplots()
ax1.scatter(psin, qpar, c=time, cmap="inferno", edgecolors="k")
ax1.set_xlim([0.98, 1.08])
ax1.set_title("#{} {}-{} ms".format(shot, tmin, tmax))
ax1.set_xlabel("psin", fontsize=14)
ax1.set_ylabel("qpar (W/cm2)", fontsize=14)
fig.tight_layout()
fig.show()
