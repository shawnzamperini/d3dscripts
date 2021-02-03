import matplotlib.pyplot as plt
import numpy as np


otf_tot = [0.220,0.274,0.249,0.297,0.112,0.165,0.060,0.063,0.268,0.034,0.028,0.032,0.012]
ps = ["A33","A34","A35","A17","A18","A19","A20","A21","A28","A22","A23","A24","A25"]
x = np.arange(0, len(otf_tot))

fig, ax = plt.subplots(figsize=(8,5))
ax.bar(x, otf_tot, color="tab:red")
ax.set_xticks(x)
ax.set_xticklabels(ps)
ax.tick_params(which='both', labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlabel("Probe Number", fontsize=16)
ax.set_ylabel("OTF Total W / Shot", fontsize=16)
fig.tight_layout()
fig.show()
