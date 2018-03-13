import atomic.atomic as atomic
import numpy as np
import pandas as pd


excel_file = "/home/shawn/d3dscripts/Data/LP_with_fit.xlsx"

# Use plunge 2 of 167195 since it gets the closest to the separatrix.
df_195_2 = pd.read_excel(excel_file,
                         sheet_name="LP Data",
                         skiprows=[0, 1],
                         names=["Time (ms)", "R (cm)", "ne (e18 m-3)", "Te (eV)",
                                "Vfl (V)", "Vp (V)", "R-Rsep (cm)"],
                         usecols=[57, 58, 59, 60, 61, 62, 63])

nes = np.abs(df_195_2["ne (e18 m-3)"].values)
Tes = np.abs(df_195_2["Te (eV)"].values)
ad = atomic.element('tungsten')
S = ad.coeffs['ionisation']
ion_coeffs = S(0, Tes, nes)
#ion_coeffs = 10**ion_coeffs

np.savetxt("ioniz_coeffs.csv", ion_coeffs, delimiter=',')
