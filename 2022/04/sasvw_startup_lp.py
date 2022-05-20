import get_lp
import matplotlib.pyplot as plt
import numpy as np


#shot = 189012
#tstart = 2000
#tend = 5500
shot = 176508
tstart = 2000
tend = 5000
lp_dict = get_lp.get_dict_of_lps(shot, tunnel=False)
data = get_lp.plot_contours(lp_dict, tstart, tend, te_vmax=75)

fig, axs = plt.subplots(2, 5)
axs = axs.flatten()
i = 0
for label in np.unique(data["labels"]):
    mask = data["labels"]==label
    time = data["time"][mask]
    te_raw = data["te_raw"][mask]
    te = data["te"][mask]

    axs[i].plot(time, te_raw, color="r", alpha=0.5)
    axs[i].plot(time, te, color="r", label=label)
    axs[i].set_title(label)

    i += 1

fig.suptitle(shot)
fig.tight_layout()
fig.show()
