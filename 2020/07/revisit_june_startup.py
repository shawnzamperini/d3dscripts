from gadata import gadata
import numpy as np
import matplotlib.pyplot as plt


def pull_avg_spred(shot, start, end, spred=None):
    """

    """
    # Load the gadata object if not supplied.
    if spred == None:
        spred = gadata("ciii_977", shot)

    # Get indicies to get average for.
    idx = np.where(np.logical_and(spred.xdata>start, spred.xdata<end))[0]
    spred_avg = spred.zdata[idx].mean()

    return spred, spred_avg

# Brief summary of the UOB parameters for each shot.
# 178346 1.2V (2-3s), 1.7V (3-4s), 2.4V (4-5s)
# 178348 2.4V (2-3s), 3.4V (3-4s), 4.4V (4-5.5s)
# 178350 1.2V (2-3s), 1.7V (3-4s), 2.4V (4-5.5s)
# 178351 2.4V (2-3s), turned off with decay in voltage.

# A time delay for signal detection in SPRED.
sig_delay = 300

# Dictionary to hold results.
startup = {"uob_346":[0, 1.2, 1.7, 2.4], "uob_348":[0, 2.4, 3.4, 4.4],
           "uob_350":[0, 1.2, 1.7, 2.4], "uob_351":[0, 2.4, 0, 0], 178346:[],
           178348:[], 178350:[], 178351:[]}

# Loop through one shot at a time.
for shot, spred_vals in startup.items():

    # The shots, which we want, are the ints in the dict.
    if type(shot) != int:
        continue

    # For these shots, just take the averages from 1900-2100 as a baseline since
    # there wasn't much of a flat-top before injection.
    spred = None
    spred, spred_bkg = pull_avg_spred(shot, 1900, 2100)
    startup[shot].append(spred_bkg)

    # Just gonna use the same times for each shot.
    for times in [(2000, 3000), (3000, 4000), (4000, 5000)]:
        spred, spred_avg = pull_avg_spred(shot, times[0]+sig_delay, times[1], spred=spred)
        startup[shot].append(spred_avg)

# Plotting commands.
ms = 14
fig, ax = plt.subplots()
ax.plot(startup["uob_346"], startup[178346], '.', ms=ms, mec='k', color="tab:purple", label=178346)
ax.plot(startup["uob_348"], startup[178348], '.', ms=ms, mec='k', color="tab:red", label=178348)
ax.plot(startup["uob_350"], startup[178350], '.', ms=ms, mec='k', color="tab:cyan", label=178350)
ax.plot(startup["uob_351"][0:2], startup[178351][0:2], '.', ms=ms, mec='k', color="tab:pink", label=178351)
ax.grid(alpha=0.5)
ax.legend(fontsize=12)
ax.set_xlabel("UOB Voltage (V)", fontsize=16)
ax.set_ylabel("SPRED 977 nm (CIII)", fontsize=16)
fig.tight_layout()
fig.show()
