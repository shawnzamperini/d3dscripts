import pretty_plots as pp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

path = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/My Slides and Sheets/2018/2018-10/avg_lams_profs.xlsx'
df = pd.read_excel(path, skiprows=[0,1], sheet_name='Peaking')
df = df[:12][['Probe', 'OTF/ITF Peaking', 'ITF/OTF Total', 'ITF/OTF Total Error']]

pnames =  df['Probe'].values
peaking = df['OTF/ITF Peaking'].values
total =   df['ITF/OTF Total'].values
total_err = df['ITF/OTF Total Error'].values

# B/C 4/5 are the reverse probes [7,9,10,11]. First make the bar chart.

for_idx = np.array([0,1,2,3,4,5,6,8])
rev_idx = np.array([7,9,10,11])
#fig = pp.pbar(y=peaking[for_idx], y2=peaking[rev_idx], bar_names=pnames[for_idx],
#              bar_names2=pnames[rev_idx])

# A good sized figure.
fig = plt.figure(figsize=(10, 7.5))
ax1 = fig.add_subplot(111)

# Remove frame lines.
ax1.spines["top"].set_visible(False)
#ax1.spines["bottom"].set_visible(False)
ax1.spines["right"].set_visible(False)
#ax1.spines["left"].set_visible(False)
ax1.set_facecolor('white')

# Axis ticks only on bottom and left.
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

# Make sure ticks are large enough to read.
ax1.tick_params(axis='both', which='both', labelsize=18)
#ax1.set_xlabel(xlabel, fontsize=fontsize, weight=weight)
ax1.set_ylabel('OTF/ITF Peaking', fontsize=26, weight='normal')

ax1.bar(for_idx, peaking[for_idx], color=pp.tableau20[6], label='Forward')
ax1.bar(rev_idx, peaking[rev_idx], color=pp.tableau20[8], label='Reverse')

ax1.legend(prop=dict(weight='normal', size=26))
ax1.set_xticks(np.arange(len(pnames)))
ax1.set_xticklabels(pnames, rotation=0)

#fig.tight_layout()
#fig.show()

fig = pp.pplot(peaking[for_idx], total[for_idx], xerr=np.full(len(peaking[for_idx]), 0.05),
               yerr=total_err[for_idx], label='Forward')
fig = pp.pplot(peaking[rev_idx], total[rev_idx], xerr=np.full(len(peaking[rev_idx]), 0.05),
               yerr=total_err[rev_idx], label='Reverse', fig=fig, color=8,
               xlabel='OTF/ITF Peaking', ylabel='ITF/OTF Total')
fig.axes[0].axhline(1.0, linestyle='--', linewidth=2, color='k', alpha=0.75)
fig.axes[0].axvline(1.0, linestyle='--', linewidth=2, color='k', alpha=0.75)
