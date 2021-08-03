from gadata import gadata
import matplotlib.pyplot as plt
import numpy as np


ga247 = gadata("RVSOUT", 167247)
ga277 = gadata("RVSOUT", 167277)


fig, ax = plt.subplots()
ax.plot(ga247.xdata, ga247.zdata, color="k", label=167247)
ax.plot(ga277.xdata, ga277.zdata, color="magenta", label=167277)
ax.set_xlim([1000, 5500])
ax.set_ylim([1.25, 1.40])
ax.fill_between(ga247.xdata, 1.32, 1.37, alpha=0.5, color="grey")
ax.set_xlabel("Time (ms)", fontsize=16)
ax.set_ylabel("Strike Point (m)", fontsize=16)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xticks([1000, 2000, 3000, 4000, 5000])
#ax.set_yticks()
ax.tick_params(labelsize=12)
ax.legend(fontsize=12)
fig.tight_layout()
fig.show()
