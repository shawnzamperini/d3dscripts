# Helper script to print out Bp and Bt along the length of an RCP plunge. In
# order to use this script you need to have pickled already the relevant EFIT
# from an OMFIT EFIT session. This can be done by:
# 1. Open OMFIT and load the EFIT module.
# 2. Enter you shot and a time during the plunge. Run EFIT with defaults,
#      choosing EFIT02 if possible.
# 3. Navigate to OMFIT['EFIT']['FILES']['gEQDSK']['AuxQuantities'] and save
#      with the following convention: shot_time_value. Ex. 184154_3000_r.
# 4. Download the four files for R, Z, Bp and Bt in d3d_work/gfile_data_for_rcp.
# 5. Change the rcp_name below and add the name of the EFIT files in chankin_const.

import pandas as pd
import numpy as np
import chankin_const
import pickle
from scipy.interpolate import interp2d
from importlib import reload


# If adding EFIT names as we go, we'll need to reload each time.
reload(chankin_const)

rcp_name = "XP184533_2"
#rcp_df = pd.read_excel("/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/rcp_master.xlsx", sheet_name=rcp_name)
rcp_df = pd.read_excel("/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/rcp_master.xlsx", sheet_name=rcp_name)

def load_pickle(path):
    with open(path, "rb") as f:
        var = pickle.load(f)
    return var

# Load in the respective EFIT data.
#root = "/mnt/c/Users/Shawn/Documents/d3d_work/gfile_data_for_rcp/"
root = "/Users/zamperini/Documents/d3d_work/gfile_data_for_rcp/"
r_efit  = load_pickle(root + chankin_const.efit_file[rcp_name] + "_r")
z_efit  = load_pickle(root + chankin_const.efit_file[rcp_name] + "_z")
bp_efit = load_pickle(root + chankin_const.efit_file[rcp_name] + "_bp")
bt_efit = load_pickle(root + chankin_const.efit_file[rcp_name] + "_bt")

f_bp = interp2d(r_efit, z_efit, bp_efit)
f_bt = interp2d(r_efit, z_efit, bt_efit)

if rcp_name[:2] == "MP":
    r = rcp_df["R (cm)"] / 100
    z = np.full(len(rcp_df), -0.188)
elif rcp_name[:2] == "XP":
    r = np.full(len(rcp_df), 1.493)
    z = rcp_df["Z (cm)"] / 100
else:
    print("Error: Did not identify probe name {}".format(rcp_name))

if rcp_name[:2] == "MP":
    bp = f_bp(r, z)[0]
    bt = f_bt(r, z)[0]
elif rcp_name[:2] == "XP":
    bp = f_bp(r, z)[:, 0]
    bt = f_bt(r, z)[:, 0]
b  = np.sqrt(bp**2 + bt**2)

print("Bp")
for i in bp:
    print(i)

print("\nBt")
for i in bt:
    print(i)
