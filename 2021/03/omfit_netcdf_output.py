import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ncpath = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/OMFITprofiles_167247_FIT.nc"
nc = netCDF4.Dataset(ncpath)

# Extract relevant data.
psin = nc.variables["psi_n"][:].data
te_ts = nc.variables["T_e"][:].mean(axis=0)
ne_ts = nc.variables["n_e"][:].mean(axis=0)
te_lp = nc.variables["T_e_LP"][:].mean(axis=0)
jsat_lp = nc.variables["n_e_LP"][:].mean(axis=0)

# Park it into a DataFrame and save as an Excel file.
df = pd.DataFrame({"psin":psin, "te_ts":te_ts, "ne_ts":ne_ts, "te_lp": te_lp,
  "jsat_lp":jsat_lp})
df.to_excel("/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/omfit_fit_ts_lp.xlsx")
