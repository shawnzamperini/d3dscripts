import pandas as pd

filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
           'Polodial_Scans/New Map Script Results/BD03_Map_Analysis.xlsx'

# Read from Excel.
df = pd.read_excel(filename, sheet_name='MapData')

# Pivot into a 2D dataFrame of total W.
df = df.pivot(index='Axial Location [mm]', columns='z Location [mm]', values='Total W')

# The mean for the average poloidal profile.
s_pol = df.mean(axis=0)

x = s_pol.index.values
y = s_pol.values / s_pol.values.max()

print("Poloidal locations:")
for i in x:
    print(i)

print("Average W:")
for i in y:
    print(i)
