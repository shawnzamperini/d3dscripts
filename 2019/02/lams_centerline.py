import pandas as pd
import numpy as np
import pretty_plots as pp
from tkinter import filedialog, Tk


# Have user pick file.
root = Tk(); root.withdraw()
print("Choose ITF File...")
xl_path_itf = filedialog.askopenfilename(filetypes=(('Excel Files', '*.xlsx'),))
print("Choose OTF File...")
xl_path_otf = filedialog.askopenfilename(filetypes=(('Excel Files', '*.xlsx'),))

xl_df_itf = pd.read_excel(xl_path_itf, sheet_name='MapData')
xl_df_otf = pd.read_excel(xl_path_otf, sheet_name='MapData')

# DataFrames to hold total counts.
total_itf = pd.DataFrame()
total_otf = pd.DataFrame()

# Pull r and z locations into their own variables.
zlocs_itf = xl_df_itf['z Location [mm]'].unique()
zlocs_otf = xl_df_otf['z Location [mm]'].unique()
rlocs = xl_df_itf['Axial Location [mm]'].unique()

# Restructure the data from the big DataFrame into the individual DataFrames.
for zloc in zlocs_itf:
    total_itf[zloc]      = xl_df_itf[xl_df_itf['z Location [mm]'] == zloc]['Total W'].values
for zloc in zlocs_otf:
    total_otf[zloc]      = xl_df_otf[xl_df_otf['z Location [mm]'] == zloc]['Total W'].values

# Set the indices to the r locations.
for df in [total_itf, total_otf]:
    df.set_index(rlocs, inplace=True)

# Get the centerline for each df.
center_itf = zlocs_itf.max() / 2.0
itf_idx = np.where(np.abs(zlocs_itf - center_itf) == np.abs(zlocs_itf - center_itf).min())[0][0]
print("ITF Centerline at {:.2f} mm. Idx = {}".format(zlocs_itf[itf_idx], itf_idx))
center_otf = zlocs_otf.max() / 2.0
otf_idx = np.where(np.abs(zlocs_otf - center_otf) == np.abs(zlocs_otf - center_otf).min())[0][0]
print("OTF Centerline at {:.2f} mm. Idx = {}".format(zlocs_otf[otf_idx], otf_idx))


y_itf = total_itf[zlocs_itf[itf_idx]].values
y_otf = total_otf[zlocs_otf[otf_idx]].values

fig = pp.pplot(rlocs, y_itf, label='ITF', fmt='-')
fig = pp.pplot(rlocs, y_otf, label='OTF', fmt='-',
               xlabel='Distance along probe (mm)',
               ylabel='LAMS counts',
               color=8, fig=fig)
