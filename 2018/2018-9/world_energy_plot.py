import pandas as pd
import pretty_plots as pp
import matplotlib.pyplot as plt



filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
           'energy_use_world.xlsx'
df = pd.read_excel(filename, skiprows=4)

fig = pp.pplot(df['Year'], df['Energy use: Total World quad Btu'], xlabel='Year',
               ylabel='World Energy Use (quad Btu)', fmt='-')

# Stack plot.

filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
           'world_energy_by_sector.xlsx'

df = pd.read_excel(filename, skiprows=4)

# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

x = df['Year']
y = [df[df.columns[1]], df[df.columns[2]], df[df.columns[3]], df[df.columns[4]], df[df.columns[5]]]

plt.rcParams['font.family'] = 'serif'
fig = plt.figure(figsize=(10, 7.5))
ax1 = fig.add_subplot(111)
ax1.stackplot(x, y, labels=['Hydro', 'Natural Gas', 'Coal', 'Nuclear', 'Other'],
              colors=tableau20[3::2])
ax1.spines["top"].set_visible(False)
#ax1.spines["bottom"].set_visible(False)
ax1.spines["right"].set_visible(False)
#ax1.spines["left"].set_visible(False)
# Axis ticks only on bottom and left.
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
# Make sure ticks are large enough to read.
ax1.tick_params(axis='both', which='both', labelsize=18)
ax1.set_xlabel('Year', fontsize=28)
ax1.set_ylabel('Energy Consumption (quad Btu)', fontsize=26)
ax1.set_xlim([2010, 2050])

handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[::-1], labels[::-1], loc='upper left', fontsize=18, framealpha=1.0)

#ax1.legend(loc='upper left', fontsize=18, framealpha=1.0)
fig.show()
