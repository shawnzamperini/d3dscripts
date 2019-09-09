import rz_to_rminrsep as omp
import numpy as np
import pandas as pd


xl_path = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2019/07/lp_with_a2.xlsx'
df = pd.read_excel(xl_path, skiprows=1)

rs = df['R (cm)'] / 100
zs = np.full(len(rs), -0.185)
times = np.arange(2500, 5000, 500)

rs_omp = omp.rz_to_rminrsep_omp(rs, zs, 167196, times)
