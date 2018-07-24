import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


filename = '/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets/2018-7/scaling.xlsx'
df = pd.read_excel(filename, sheet_name='Scaling Law 1', skiprows=12)

# Drop everything after, except the three unique ones.
df.drop([0, 1, 2, 3, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24], inplace=True)

ratios = np.array(df['Ratio'].values, dtype=np.float64)
ratio_err = np.array(df['Error'].values, dtype=np.float64)
lambs = np.array(df['Density Fall Off (cm)'].values, dtype=np.float64)
lambs_err = np.array(df['Error (cm)'].values, dtype=np.float64)

plt.style.use('seaborn')
font = {'fontsize':24}
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.errorbar(lambs, ratios, ratio_err, lambs_err, 'k.', markeredgewidth=1, capsize=5, ms=15)
ax1.set_xlabel(r'$\mathrm{\lambda_{ne}}$ (cm)',  font)
ax1.set_ylabel('ITF/OTF Total W Ratio', font)
ax1.tick_params(labelsize=22)
fig.tight_layout()
plt.show()
