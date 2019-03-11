import pandas as pd
import numpy as np
import pretty_plots as pp
from tkinter import filedialog, Tk


# Have user pick file.
root = Tk(); root.withdraw()
xl_path = filedialog.askopenfilename(filetypes=(('Excel Files', '*.xlsx'),))
print(xl_path)

# Load in all the data into a DataFrame.
xl_df = pd.read_excel(xl_path, sheet_name='MapData')

# DataFrame to hold the rearranged data.
total      = pd.DataFrame()
shelf_frac = pd.DataFrame()
floor_frac = pd.DataFrame()
shelf_tot  = pd.DataFrame()
floor_tot  = pd.DataFrame()

# Pull r and z locations into their own variables.
zlocs = xl_df['z Location [mm]'].unique()
rlocs = xl_df['Axial Location [mm]'].unique()

# Restructure the data from the big DataFrame into the individual DataFrames.
for zloc in zlocs:
    total[zloc]      = xl_df[xl_df['z Location [mm]'] == zloc]['Total W'].values
    shelf_frac[zloc] = xl_df[xl_df['z Location [mm]'] == zloc]['Shelf Fraction'].values
    floor_frac[zloc] = xl_df[xl_df['z Location [mm]'] == zloc]['Floor Fraction'].values
    shelf_tot[zloc]  = xl_df[xl_df['z Location [mm]'] == zloc]['Shelf Total'].values
    floor_tot[zloc]  = xl_df[xl_df['z Location [mm]'] == zloc]['Floor Total'].values

# Set the indices to the r locations.
for df in [total, shelf_frac, floor_frac, shelf_tot, floor_tot]:
    df.set_index(rlocs, inplace=True)

# Function for plotting a DataFrame, including option to clip to a top value.
def plot_contourf(df, title, clip_value=None):
    x = df.index.values
    y = df.columns.values
    X, Y = np.meshgrid(x, y)
    Z = df.values.T

    if clip_value:
        Z = np.clip(Z, 0, clip_value)
    else:
        Z = np.clip(Z, 0, 1e13)

    fig = pp.ppcontourf(X, Y, Z, xlabel='Radial (mm)', ylabel='Poloidal (mm)',
                        cbarlabel=title, extend='max')

# Inital plots.
plot_contourf(total, 'Total W')
plot_contourf(shelf_frac, 'Shelf Fraction', 1.0)
plot_contourf(floor_frac, 'Floor Fraction', 1.0)

prompt = 'Choose range to get average fraction (separated by commas): '
min_val, max_val = [float(x) for x in input(prompt).split(',')]

df1 = floor_frac.loc[floor_frac.index > min_val]
df2 = df1.loc[df1.index < max_val]

print('Average floor fraction = {:.2f}'.format(df2.mean().mean()))
