import get_ts as ts
from scipy.interpolate import Rbf
from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt

shot = 167277
tmin = 2500
tmax = 4000

# Get core TS data (i.e. the dat at the top of the device).
core = ts.run_script(shot, 'core', tmin=2500, tmax=4500, tstep=300)

# Get the divertor TS data (i.e. close to the target).
div = ts.run_script(shot, 'divertor', tmin=2500, tmax=4500, tstep=300)

# Get average psin of divertor TS closest to target.
div_psin = div['psins']['avg_psins'][0]

# Get the raw TS data for the divertor TS.
div_raw_time = div['temp']['X']
div_raw_temp = div['temp']['Y'][0]

# Filter out the spikes from ELMs by just removing data above a threshold.
thresh = 40
div_noelm_time = div_raw_time[np.where(div_raw_temp<thresh)]
div_noelm_temp = div_raw_temp[np.where(div_raw_temp<thresh)]

# Filter the data using a savgol filter.
div_filt_temp = savgol_filter(div_noelm_temp, 51, 3)

# Finally get the average temp in the requested time range.
div_avg_temp = np.mean(div_filt_temp[np.where(np.logical_and(div_noelm_time>tmin, div_noelm_time<tmax))])

# Now we need the average temperature upstream from the core TS at the same psin.
core_psins  = core['psins']['avg_psins']
closest_idx = np.abs(core_psins-div_psin).argmin()
core_psin   = core_psins[closest_idx]

# Get the average temperature at this psin.
core_raw_time = core['temp']['X']
core_raw_temp = core['temp']['Y'][closest_idx]

# Same cleaning up process as for the divertor.
thresh = 110
core_noelm_time = core_raw_time[np.where(core_raw_temp<thresh)]
core_noelm_temp = core_raw_temp[np.where(core_raw_temp<thresh)]

core_filt_temp = savgol_filter(core_noelm_temp, 51, 3)

core_avg_temp = np.mean(core_filt_temp[np.where(np.logical_and(core_noelm_time>tmin, core_noelm_time<tmax))])

print("Based on chord nearest target:")
print("  Divertor psin used:           {:.4f}".format(div_psin))
print("  Closest psin from core:       {:.4f}".format(core_psin))
print("  Divertor average temperature: {:.4f}".format(div_avg_temp))
print("  Core average temperature:     {:.4f}".format(core_avg_temp))
print("  Upstream/downstream ratio:    {:.4f}".format(core_avg_temp/div_avg_temp))

# Try the values closest to the separatrix.
div_chord_near_sep  = np.max(np.where(div['psins']['avg_psins']>1.0))
core_chord_near_sep = np.max(np.where(core['psins']['avg_psins']>1.0))

div_raw_time = div['temp']['X']
div_raw_temp = div['temp']['Y'][div_chord_near_sep]
div_thresh = 100
div_noelm_time = div_raw_time[np.where(div_raw_temp<div_thresh)]
div_noelm_temp = div_raw_temp[np.where(div_raw_temp<div_thresh)]
div_filt_temp = savgol_filter(div_noelm_temp, 51, 3)
div_avg_temp = np.mean(div_filt_temp[np.where(np.logical_and(div_noelm_time>tmin, div_noelm_time<tmax))])

core_psin   = core_psins[core_chord_near_sep]
core_raw_time = core['temp']['X']
core_raw_temp = core['temp']['Y'][core_chord_near_sep]
core_thresh = 110
core_noelm_time = core_raw_time[np.where(core_raw_temp<core_thresh)]
core_noelm_temp = core_raw_temp[np.where(core_raw_temp<core_thresh)]
core_filt_temp = savgol_filter(core_noelm_temp, 51, 3)
core_avg_temp = np.mean(core_filt_temp[np.where(np.logical_and(core_noelm_time>tmin, core_noelm_time<tmax))])

print("Based on chord nearest separatrix:")
print("  Divertor psin used:           {:.4f}".format(div_psin))
print("  Closest psin from core:       {:.4f}".format(core_psin))
print("  Divertor average temperature: {:.4f}".format(div_avg_temp))
print("  Core average temperature:     {:.4f}".format(core_avg_temp))
print("  Upstream/downstream ratio:    {:.4f}".format(core_avg_temp/div_avg_temp))

# Try getting the upstream/downstream ratio for each divertor chord location.
print("Div  Te     Core Te")
for chord in range(0, 8):
    div_psin = div['psins']['avg_psins'][chord]
    div_raw_time = div['temp']['X']
    div_raw_temp = div['temp']['Y'][chord]
    div_thresh = 100
    div_noelm_time = div_raw_time[np.where(div_raw_temp<div_thresh)]
    div_noelm_temp = div_raw_temp[np.where(div_raw_temp<div_thresh)]
    div_filt_temp = savgol_filter(div_noelm_temp, 51, 3)
    div_avg_temp = np.mean(div_filt_temp[np.where(np.logical_and(div_noelm_time>tmin, div_noelm_time<tmax))])

    closest_idx = np.abs(core_psins-div_psin).argmin()
    core_psin   = core_psins[closest_idx]
    core_raw_time = core['temp']['X']
    core_raw_temp = core['temp']['Y'][closest_idx]
    core_thresh = 110
    core_noelm_time = core_raw_time[np.where(core_raw_temp<core_thresh)]
    core_noelm_temp = core_raw_temp[np.where(core_raw_temp<core_thresh)]
    core_filt_temp = savgol_filter(core_noelm_temp, 51, 3)
    core_avg_temp = np.mean(core_filt_temp[np.where(np.logical_and(core_noelm_time>tmin, core_noelm_time<tmax))])
    print("{:4} {:5.3f} {:4} {:5.3f}".format(chord, div_avg_temp, closest_idx, core_avg_temp))

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(div_raw_time, div_raw_temp, 'k', label='Raw')
ax1.plot(div_noelm_time, div_filt_temp, 'r', label='Filtered')
ax1.legend()
ax1.set_ylabel('Div Te (eV)')
ax1.set_xlim([0, 6000])
ax1.set_ylim([0, 500])
ax1.axhline(y=div_thresh, linestyle='--')
ax2.plot(core_raw_time, core_raw_temp, 'k', label='Raw')
ax2.plot(core_noelm_time, core_filt_temp, 'r', label='Filtered')
ax2.legend()
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('Core Te (eV)')
ax2.set_xlim([0, 6000])
ax2.set_ylim([0, 500])
ax2.axhline(y=core_thresh, linestyle='--')
fig.tight_layout()
fig.show()
