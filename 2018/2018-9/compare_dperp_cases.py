import matplotlib.pyplot  as plt
import pandas             as pd
import numpy              as np


# Location of Excel file with "TOTAL DEPOSITION" table copy/pasted into it.
filename7 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/test24_totdep.xlsx'
filename8 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/test25_totdep.xlsx'

norm = True

# Load the data into a Dataframe, drop the last row that has junk in it, and
# then rename the column that is a string '0.1' to a float 0.0 (not sure why
# this happens).
print("Loading file...")
df7 = pd.read_excel(filename7, index_col=0)
df7.drop(df7.columns[-1], axis=1, inplace=True)
df7.rename({'0.1':0.0, '0.2':0.0}, axis=1, inplace=True)
df8 = pd.read_excel(filename8, index_col=0)
df8.drop(df8.columns[-1], axis=1, inplace=True)
df8.rename({'0.1':0.0, '0.2':0.0}, axis=1, inplace=True)

# Get the data at P=0 (so it will be radial centerline plots).
x7r = df7[0].index.values
y7r = df7[0].values * -1
x8r = df8[0].index.values
y8r = df8[0].values * -1

# Now get the data at X = 0 (so poloidal plots). Will use 0.00029925 as X=0.
r_val = 0.0032918
x7p = df7.loc[r_val].index.values
y7p = df7.loc[r_val].values * -1
x8p = df8.loc[r_val].index.values
y8p = df8.loc[r_val].values * -1

if norm:
    y7r = y7r/y7r.max()
    y8r = y8r/y8r.max()
    y7p = y7p/y7p.max()
    y8p = y8p/y8p.max()


# Plotting commands.
plt.style.use('seaborn')
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x7r, y7r, 'C2.', label='Test 24')
ax1.plot(x8r, y8r, 'C3.', label='Test 25')
ax1.set_xlabel('x', fontsize=24)
ax1.set_ylabel('Deposition (arbitrary units)', fontsize=24)
ax1.legend(fontsize=24)
fig.tight_layout()
fig.show()

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x7p, y7p, 'C2', lw=5, label='Test 24')
ax1.plot(x8p, y8p, 'C3', lw=5, label='Test 25')
ax1.set_xlabel('Poloidal (m)', fontsize=24)
ax1.set_ylabel('Deposition (arbitrary units)', fontsize=24)
ax1.set_xlim([-0.025, 0.025])
ax1.legend(fontsize=24)
fig.tight_layout()
fig.show()
