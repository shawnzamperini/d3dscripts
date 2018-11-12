import gadata
import MDSplus
import numpy as np
import pretty_plots as pp
import matplotlib.pyplot as plt


shot = int(input('Shot: '))

conn = MDSplus.Connection('localhost')
fs04 = gadata.gadata('FS04', shot, connection=conn)

dz = np.zeros(fs04.xdata.shape, np.float)
dz[0:-1] = np.diff(fs04.zdata) / np.diff(fs04.xdata)
dz[-1] = (fs04.zdata[-1] - fs04.zdata[-2]) / (fs04.xdata[-1] - fs04.xdata[-2])

#fig = pp.pplot(fs04.xdata, fs04.zdata, '-', lw=1, xlabel='Time (ms)', ylabel='FS04')
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(fs04.xdata, fs04.zdata)
ax2.plot(fs04.xdata, dz)
fig.show()

min_time = float(input('Min. Time: '))
max_time = float(input('Max. Time: '))
min_time_idx = np.where(fs04.xdata > min_time)[0].min()
max_time_idx = np.where(fs04.xdata < max_time)[0].max()
thresh = 0.75e17
mask = np.array(dz > thresh)

# Set all values behind min_time and ahead of max_time in mask to False.
mask[:min_time_idx] = False
mask[max_time_idx:] = False

# Prevent double counting of ELM by allowing a grace period of elm_length before
# considering another ELM.
elm_length = 2
time_res = (fs04.xdata[-1] - fs04.xdata[0]) / fs04.xdata.size
idx_per_elm = int(elm_length/time_res)
count = 0
for mask_val in mask:
    if mask_val:
        mask[count+1:count+idx_per_elm] = False
    count += 1

# Find the ELM frequency of the selected time range.
n_elms = np.sum(mask)
elm_freq = n_elms / (max_time - min_time)
print("ELM Frequency = {:.3f} Hz".format(elm_freq*1000))

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(fs04.xdata, fs04.zdata)
#ax1.set_xlim([2500, 2600])
ax2.plot(fs04.xdata, dz)
#ax2.set_xlim([2500, 2600])
ax2.axhline(thresh, color='r', lw=2)
ax2.set_ylim([-3e17, 3e17])
for line in fs04.xdata[mask]:
    ax1.axvline(line, linestyle='--', color='r')
    ax2.axvline(line, linestyle='--', color='r')
ax2.set_xlabel('Time (ms)')
ax1.set_ylabel('FS04 (ph/cm2/sr/s)')
ax2.set_ylabel('Derivative of FS04')
fig.show()

print('Plot time range (q to quit):')
while True:
    try:
        new_min = float(input('  Time min:'))
        new_max = float(input('  Time max:'))
        ax1.set_xlim([new_min, new_max])
        ax2.set_xlim([new_min, new_max])
        fig.show()
    except:
        break
