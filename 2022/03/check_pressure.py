import get_lp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


shot = 184267

keep_names = np.array(["probe {}".format(n) for n in np.arange(75, 102)])

lp_dict = get_lp.plot_lps(shot, 2700, 3700, "psin", 10, tunnel=False, up=True)
keep = np.full(len(lp_dict["psin"]), False)
for i in range(0, len(keep)):
    if lp_dict["pnames"][i] in keep_names:
        keep[i] = True
psin = np.array(lp_dict["psin"])[keep]
ne = np.array(lp_dict["ne (cm-3)"])[keep] * 1e6
te = np.array(lp_dict["Te (eV)"])[keep]
mb = 2.0 * 931.49 * 10**6 / ((3*10**8)**2)   # eV * s2 / m2
cs = np.sqrt((te + te) / mb)
mb_kg = 3.32e-27
ptarg = ne * te * 1.619e-19 + ne * mb_kg * np.square(cs)   # Pascals

rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"
rcp = pd.ExcelFile(rcp_path)
mp1 = rcp.parse("XP{}_1".format(shot))
pomp = mp1["ne (1e18 m-3)"] * 1e18 * mp1["Te (eV)"] * 1.619e-19 + mp1["ne (1e18 m-3)"] * 1e18 * mb_kg * np.square(mp1["Vflow (m/s)"])

fig, ax = plt.subplots()
ax.plot(psin, ptarg, label="LP", color="tab:red", lw=3)
ax.plot(mp1["Psin"], pomp, label="XRCP", color="tab:purple", lw=3)
ax.set_xlabel("Psin")
ax.set_ylabel("Pressure (Pa)")
ax.grid()
ax.legend()
fig.tight_layout()
fig.show()
