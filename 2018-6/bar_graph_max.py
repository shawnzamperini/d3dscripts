import matplotlib.pyplot as plt
import pandas as pd
import numpy  as np


filename = '/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets/2018 - 5/max_w_vs_flux.xlsx'
df = pd.read_excel(filename, sheet_name='Bar Graph', usecols='A:C')

ind = np.arange(len(df['Max*100']))

font = {'fontsize' : 24,
        'weight'   : 'bold'}
plt.style.use('seaborn')
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.barh(ind, df['Max*100'], color='black')
ax1.get_yaxis().set_visible(False)
ax1.set_xlabel(r'$\mathrm{\bf{Max\ W\ Content\ (10^{12}\ cm^{-2})}}}$', font)
ax1.set_xscale('log')
ax1.set_xlim([1, 200])
ax1.grid(True, which='minor', axis='x')
ax1.tick_params(axis='x', which='major', labelsize=20)
fig.tight_layout()
fig.show()
