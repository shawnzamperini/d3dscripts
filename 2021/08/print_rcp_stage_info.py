# This script uses the gfile info to calculate R-Rsep along the length of the
# plunges, as well as print out the Bp and Bt values which are all needed to
# calculate Chankin's formula. The below script can be copied into OMFIT after
# running the EFIT module for the specific case to save the gfile info into
# my gfile folder.
"""
import pickle

g_root = OMFIT["EFIT"]["FILES"]["gEQDSK"]
aux_root = g_root["AuxQuantities"]
shot = g_root["CASE"][3][2:]
time = g_root["CASE"][4][:-2]

# Wrap everything into one dictionary.
gfile = {}
gfile["R"] = aux_root["R"]
gfile["Z"] = aux_root["Z"]
gfile["Bp"] = aux_root["Bp"]
gfile["Bt"] = aux_root["Bt"]
gfile["RBBBS"] = g_root["RBBBS"]
gfile["ZBBBS"] = g_root["ZBBBS"]

# Deploy as pickle funniest shit I ever saw.
fname = "/home/zamperinis/gfiles/" + shot + "_" + str(int(time))
with open(fname, "wb") as f:
    pickle.dump(gfile, f)
"""

import pickle
import numpy  as np
import pandas as pd
from scipy.interpolate import interp1d, interp2d

rcp_name = "MP184179_1"

rcp_efits = {"XP184533_2":"184533_3340", "MP184179_1":"184179_1600"}


# Load RCP data to pull out the locations.
rcp_df = pd.read_excel("/Users/zamperini/Google Drive/My Drive/Research/" + \
  "Data/rcp_data/rcp_master.xlsx", sheet_name=rcp_name)

# Load the pickled dictionary.
gfile_path = "/Users/zamperini/Documents/d3d_work/gfile_data_for_rcp/" + \
  "{}".format(rcp_efits[rcp_name])
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)

# Poloidal and toroidal fields.
f_bp = interp2d(gfile["R"], gfile["Z"], gfile["Bp"])
f_bt = interp2d(gfile["R"], gfile["Z"], gfile["Bt"])

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

# Calculate distance from separatrix. Limit to just portions of the LCFS
# according to the RCP since it isn't a well-defined problem (fails the straight
# line test).
if rcp_name[:2] == "MP":
    mask = gfile["RBBBS"] > 2.00
    f_lcfs_r = interp1d(gfile["ZBBBS"][mask], gfile["RBBBS"][mask])
    rsep = f_lcfs_r(z[0])
    dist = r - rsep
elif rcp_name[:2] == "XP":
    mask = gfile["ZBBBS"] < -0.80
    f_lcfs_z = interp1d(gfile["RBBBS"][mask], gfile["ZBBBS"][mask])
    zsep = f_lcfs_z(r[0])
    dist = zsep - z  # Could do abs, but do this so it's positive when we want it to be.

# Print out everything we just did for copy/paste.
print("Bp")
for i in bp:
    print(i)

print("\nBt")
for i in bt:
    print(i)

print("\nDistance from separatrix (m)")
for i in dist:
    print(i)
