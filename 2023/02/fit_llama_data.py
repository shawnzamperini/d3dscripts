import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt, savgol_filter
from scipy.optimize import curve_fit



def exp_fit(x, a, b, c):
    return a * np.exp(b * x) + c

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9, 4), sharex=True, sharey=True)
ax1.set_yscale("log")
ax1.set_xlim([0.98, 1.11])
ax1.set_ylim([1e14, 1e20])
for shot in [190484, 190485, 190486]:

    if shot == 190484:
        ax = ax1
    elif shot == 190485:
        ax = ax2
    elif shot == 190486:
        ax = ax3

    # LLAMA data.
    llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_{}_.npz".format(shot)
    llama = np.load(llama_path)
    lpsin = llama["psi_n"]
    lrho = np.square(lpsin)  # According to Florian L.
    lneut = llama["nDens_LFS"]
    lneut_err = llama["nDens_LFS_err"]

    # Filter out spikes.
    lneut_filt = medfilt(lneut, 51)
    # neut_filt = savgol_filter(neut, 51, 3)

    # Fit to exponential that increases outwards. Restrict to the SOL.
    mask = np.logical_and(lpsin > 0.98, lpsin < 1.05)
    psin = lpsin[mask]
    neut = lneut_filt[mask]
    popt, pcov = curve_fit(exp_fit, psin - 1, neut / 1e16, maxfev=5000, p0=(1, 5, 1))
    neut_fit = exp_fit(psin - 1, *popt) * 1e16

    # Print out fit.
    print("{:}: {:.2e} exp({:.2e}*x) + {:.2e}".format(shot, popt[0] * 1e16, popt[1], popt[2] * 1e16))

    ax.axvline(1.0, color="k")
    ax.plot(lpsin, lneut, alpha=0.2, color="tab:red")
    ax.plot(lpsin, lneut_filt, color="tab:red", lw=3)
    ax.plot(psin, neut_fit, color="k", lw=3)
    ax.grid()

fig.tight_layout()
fig.show()