import get_lp
import matplotlib.pyplot as plt
import sys
from scipy.signal import medfilt
from scipy.interpolate import interp1d
from gadata import gadata
import MDSplus as mds
import numpy as np


conn = mds.Connection("atlas.gat.com")

# Values of when I-mode occurred according to Amanda.
imode = {189384: [2700, 4000],
         189376: [[3400, 3750], [2800, 3100]],
         189377: [2800, 3200],
         189379: [2800, 3200],
         189380: [[3400, 3700], [2600, 3000]],
         189381: [2700, 3000],
         189382: [3100, 4300],
         189383: [3100, 3300]}

def run(shot, pnum, tmin, tmax, window=23):
    lpdict = get_lp.get_dict_of_lps(shot, tunnel=False)
    pname = "probe {}".format(pnum)
    time = lpdict[pname]["time"]
    temp = lpdict[pname]["temp"]
    heat = lpdict[pname]["heatflux"]

    # Only keep data within designated time range.
    mask = np.logical_and(time >= tmin, time <= tmax)
    time = time[mask]
    temp = temp[mask]
    heat = heat[mask]

    # The stored energy, within time range. We pad a little before/after so that we can ensure we span the entire
    # interpolation range needed below.
    ga_obj = gadata("WMHD", shot, connection=conn)
    wmhd_t = ga_obj.xdata
    wmhd = ga_obj.zdata
    mask = np.logical_and(wmhd_t >= tmin - 25, wmhd_t <= tmax + 25)
    wmhd_t = wmhd_t[mask]
    wmhd = wmhd[mask]

    # The injected power and total core radiated power.
    ga_obj = gadata("PINJ", shot, connection=conn)
    pinj_t = ga_obj.xdata
    pinj = ga_obj.zdata * 1e-3  # kW to MW
    ga_obj = gadata("PRAD_CORE", shot, connection=conn)
    prad_t = ga_obj.xdata
    prad = ga_obj.zdata * 1e-6  # W to MW

    # Median filter with a good window.
    heat_med = medfilt(heat, window)
    temp_med = medfilt(temp, window)

    # Interpolate onto a common X-axis.
    f_wmhd = interp1d(wmhd_t, wmhd)
    f_prad = interp1d(prad_t, prad)
    f_pinj = interp1d(pinj_t, pinj)
    wmhd_int = f_wmhd(time)
    psol = f_pinj(time) - f_prad(time)
    psol_med = medfilt(psol, window)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    if len(np.shape(imode[shot])) > 1:
        for lims in imode[shot]:
            ax1.axvspan(lims[0], lims[1], color="tab:red", alpha=0.3)
            ax2.axvspan(lims[0], lims[1], color="tab:red", alpha=0.3)
    else:
        ax1.axvspan(imode[shot][0], imode[shot][1], color="tab:red", alpha=0.3)
        ax2.axvspan(imode[shot][0], imode[shot][1], color="tab:red", alpha=0.3)

    ax1.plot(time, heat, color="k")
    ax1.plot(time, heat_med, color="r")
    ax1.set_ylabel("Heat Flux (W/cm2)")
    ax1.set_ylim(0, 2)

    ax2.plot(time, psol_med, color="k")
    ax2.set_ylabel("PSOL (MW)")

    ax2.set_xlabel("Time (ms)")
    ax1.set_title(shot)
    fig.tight_layout()
    fig.show()

    # Separate out the during and not during I-mode data.
    imode_mask = np.full(len(psol_med), False)
    if len(np.shape(imode[shot])) > 1:
        for lims in imode[shot]:
            imode_mask[np.logical_and(time >= lims[0], time <= lims[1])] = True
    else:
        imode_mask[np.logical_and(time >= imode[shot][0], time <= imode[shot][1])] = True

    fig, ax1 = plt.subplots(figsize=(5, 4))
    ax1.scatter(psol_med[~imode_mask], heat_med[~imode_mask], s=30, color="tab:red", label="L-Mode", edgecolor="k", marker="^")
    ax1.scatter(psol_med[imode_mask], heat_med[imode_mask], s=30, color="tab:cyan", label="I-Mode", edgecolor="k", marker="^")
    ax1.set_xlabel("PSOL (MW)")
    ax1.set_ylabel("{} Heat Flux (W/cm2)".format(pnum))
    ax1.legend()
    ax1.grid(alpha=0.3)
    fig.tight_layout()
    fig.show()


# Probe 33 was on the strike point for Amanda's I-mode experiments.
if __name__ == "__main__":
    shot = int(sys.argv[1])
    pnum = int(sys.argv[2])
    tmin = float(sys.argv[3])
    tmax = float(sys.argv[4])
    run(shot, pnum, tmin, tmax)
