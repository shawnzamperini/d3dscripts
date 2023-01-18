# A hopefully simple script to create a periodic spike with an exponential decay.
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

n_per_spike = 100

path = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190440_1.tab"
rcp = pd.read_csv(path, delimiter="\t")
ne_meas = []
avg_nes = []
rs = []
vr_meas = []
avg_vrs = []
square_vrs = []
meas_vrs = rcp["Vr(m/s)"]
for i in range(0, len(rcp)):
    rs.append(rcp["R(cm)"].iloc[i])
    f = int(rcp["Npeaks"].iloc[i] / 0.005)
    print(f)
    t = np.linspace(0, 1, f*n_per_spike)
    ne_meas.append(rcp["Ne (e18m-3)"].iloc[i] * 1e18)
    ne = np.zeros(len(t))
    vr = np.zeros(len(t))
    ne0 = 1e19
    vr0 = rcp["Vr(m/s)"].iloc[i]
    #decay = rcp["D_rad(cm)"].iloc[i] / 100 / rcp["Vr(m/s)"].iloc[i]
    decay = rcp["T_blob(e-6s)"].iloc[i] * 1e-6

    for i in range(0, f):
        tbase = np.linspace(0, 1/f, n_per_spike)
        ne[i*n_per_spike:(i+1)*n_per_spike] = ne0 * np.exp(-tbase / decay)
        vr[i * n_per_spike:(i + 1) * n_per_spike] = vr0 * np.exp(-tbase / decay)

    square_vrs.append(vr0 * decay * f)

    avg_ne = np.trapz(ne, t)
    avg_vr = np.trapz(vr, t)
    print("{:.2e}".format(avg_vr))
    avg_nes.append(avg_ne)
    avg_vrs.append(avg_vr)
print("Avg. avg. vr: {:.2f}".format(np.mean(avg_vrs)))

fig, ax = plt.subplots()
ax.plot(rs, meas_vrs)
ax.plot(rs, avg_vrs)
ax.plot(rs, square_vrs)
fig.tight_layout()
fig.show()