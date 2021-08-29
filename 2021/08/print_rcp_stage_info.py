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

rcp_name = "XP184154_2"

rcp_efits = {"XP184533_2":"184533_3340", "MP184179_1":"184179_1600",
             "XP184266_2":"184266_3340", "XP184528_1":"184528_1580",
             "XP184264_1":"184264_1530", "XP184531_1":"184531_1530",
             "XP184533_1":"184533_1530", "MP184179_2":"184179_3500",
             "XP184266_1":"184266_1530", "XP184528_2":"184528_3380",
             "XP184264_2":"184264_2940", "XP184531_2":"184531_3340",
             "XP184183_1":"184183_1540", "XP184535_2":"184535_3340",
             "MP184182_1":"184182_1620", "MP184530_2":"184530_3400",
             "MP187108_1":"187108_1620", "XP184178_2":"184178_3450",
             "MP184532_1":"184532_1610", "XP184537_1":"184537_1540",
             "MP187111_1":"187111_1620", "MP184267_1":"184267_1610",
             "MP184529_2":"184529_3400", "MP184530_1":"184530_1610",
             "XP184535_1":"184535_1540", "MP184182_2":"184182_3500",
             "MP187108_2":"187108_3410", "XP184178_1":"187178_1540",
             "XP184537_2":"184537_3340", "XP184182_1":"184182_1620",
             "MP184527_1":"184527_1600", "XP184527_1":"184527_1600",
             "XP184154_1":"184154_2940", "XP184154_2":"184154_4430"}


# Load RCP data to pull out the locations.
rcp_df = pd.read_excel("/Users/zamperini/Google Drive/My Drive/Research/" + \
  "Data/rcp_data/rcp_master.xlsx", sheet_name=rcp_name)

# Load the pickled dictionary.
#gfile_path = "/Users/zamperini/Documents/d3d_work/gfile_data_for_rcp/" + \
#  "{}".format(rcp_efits[rcp_name])
gfile_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/gfile_data_for_rcp/" + \
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
