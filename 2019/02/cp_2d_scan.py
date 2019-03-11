import pandas as pd
import numpy as np
import pretty_plots as pp
from tkinter import filedialog, Tk


root = Tk(); root.withdraw()
xl_path = filedialog.askopenfilename(filetypes=(('Excel Files', '*.xlsx'),))

print(xl_path)

df_total  = pd.read_excel(xl_path, sheet_name='Total', index_col=0)
try:
    df_ef182    = pd.read_excel(xl_path, sheet_name='EF182', index_col=0)
except:
    print("No sheet named EF182")
try:
    df_shelf_frac = pd.read_excel(xl_path, sheet_name='Shelf Fraction', index_col=0)
except:
    print("No sheet named Shelf Fraction")
try:
    df_shelf_tot  = pd.read_excel(xl_path, sheet_name='Shelf Total', index_col=0)
except:
    print("No sheet named Shelf Total")
try:
    df_floor_tot = pd.read_excel(xl_path, sheet_name='Floor Total', index_col=0)
except:
    print("No sheet named Floor Total")
try:
    df_w182 = pd.read_excel(xl_path, sheet_name='W182', index_col=0)
except:
    print("No sheet named W182")
try:
    df_w184 = pd.read_excel(xl_path, sheet_name='W184', index_col=0)
except:
    print("No sheet named W184")


def plot_contourf(df, title, clip_value=None):
    x = df.index.values
    y = df.columns.values
    X, Y = np.meshgrid(x, y)
    Z = df.values.T

    if clip_value:
        Z = np.clip(Z, 0, clip_value)

    fig = pp.ppcontourf(X, Y, Z, xlabel='Radial (mm)', ylabel='Poloidal (mm)',
                        cbarlabel=title)

plot_contourf(df_total, 'Total')
#plot_contourf(df_ef182, 'EF182')
plot_contourf(df_shelf_frac, 'Shelf Fraction', clip_value=1.0)
plot_contourf(1-df_shelf_frac, 'Floor Fraction', clip_value=1.0)
plot_contourf(df_shelf_tot, 'Shelf Total')
plot_contourf(df_floor_tot, 'Floor Total')
#plot_contourf(df_w182, 'W182')
#plot_contourf(df_w184, 'W184')

while True:
    try:
        clip_value = float(input("Clip value for Total: "))
        plot_contourf(df_total, 'Total', clip_value)
        plot_contourf(df_shelf_tot, 'Shelf Total', clip_value)
        plot_contourf(df_floor_tot, 'Floor Total', clip_value)
        #plot_contourf(df_w182, 'W182', clip_value)
        #plot_contourf(df_ef182, 'EF182', clip_value)
        #plot_contourf(df_w184, 'W184', clip_value)
    except:
        break
