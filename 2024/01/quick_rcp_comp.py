import pandas as pd
import matplotlib.pyplot as plt


rcp1 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP190412_2.tab", delimiter="\t")
rcp2 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP190440_2.tab", delimiter="\t")
rcp_r1 = rcp1["R(cm)"] / 100
rcp_r2 = rcp2["R(cm)"] / 100
rcp_ne1 = rcp1["Ne(E18 m-3)"] * 1e18
rcp_ne2 = rcp2["Ne(E18 m-3)"] * 1e18
rcp_te1 = rcp1["Te(eV)"]
rcp_te2 = rcp2["Te(eV)"]
rcp_m1 = rcp1["Machn"]
rcp_m2 = rcp2["Machn"]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(5, 8))

ax1.plot(rcp_r1, rcp_ne1, label="190412")
ax1.plot(rcp_r2, rcp_ne2, label="190440")
ax2.plot(rcp_r1, rcp_te1)
ax2.plot(rcp_r2, rcp_te2)
ax3.plot(rcp_r1, rcp_m1)
ax3.plot(rcp_r2, rcp_m2)
ax1.legend()

fig.tight_layout()
fig.show()