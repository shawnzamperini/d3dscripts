import matplotlib.pyplot as plt
import numpy as np


r = np.linspace(-3, 3, 100)

def get_f(x, birth_loc, lamb, max_freq=1000):

    f = max_freq * np.exp((x-birth_loc)/lamb)
    f[f > max_freq] = max_freq
    return f


fblob1 = get_f(r, -1, 1)
fhole1 = get_f(r, -1, -1)
fblob2 = get_f(r, 0, 1)
fhole2 = get_f(r, 0, -1)
fblob3 = get_f(r, 1, 1)
fhole3 = get_f(r, 1, -1)

lw = 3
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9, 3), sharex=True, sharey=True)

for ax in [ax1, ax2, ax3]:
    ax.axvline(0, color="k", linestyle="-")
    ax.set_xlabel("R-Rsep @ OMP (cm)")
ax1.set_ylabel("Frequency (Hz)")

ax1.axvline(-1, color="k", linestyle="--")
ax1.plot(r, fblob1, label="Blobs", color="tab:red", lw=lw)
ax1.plot(r, fhole1, label="Holes", color="tab:purple", lw=lw)
ax1.set_title("Less core contamination")

ax2.axvline(0, color="k", linestyle="--")
ax2.plot(r, fblob2, label="Blobs", color="tab:red", lw=lw)
ax2.plot(r, fhole2, label="Holes", color="tab:purple", lw=lw)
ax2.set_title("----->")

ax3.axvline(1, color="k", linestyle="--")
ax3.plot(r, fblob3, label="Blobs", color="tab:red", lw=lw)
ax3.plot(r, fhole3, label="Holes", color="tab:purple", lw=lw)
ax3.set_title("More core contamination")

ax1.legend()


fig.tight_layout()
fig.show()