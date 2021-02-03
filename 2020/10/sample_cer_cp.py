import matplotlib.pyplot as plt
import numpy as np


mu = 1
sd = 0.2

cp_x = np.linspace(1.5, 5, 30)
cp_y = np.exp(-cp_x)
cp_noise = np.random.normal(mu, sd, 30)
cp_y *= cp_noise

cer_x = np.linspace(0, 2, 7)
cer_y = np.exp(-cer_x)
cer_noise = np.random.normal(mu, sd, 7)
cer_y *= cer_noise

fig, ax = plt.subplots()
ax.plot(cp_x, cp_y, color="tab:purple", lw=3)
ax.errorbar(x=cer_x, y=cer_y, xerr=0.2, yerr=0.1, ms=8, fmt='.', color="tab:red")
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel("Distance from Separatrix", fontsize=16)
ax.set_ylabel("Carbon Density", fontsize=16)
fig.tight_layout()
fig.show()
