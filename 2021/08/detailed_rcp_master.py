import pandas as pd
import numpy  as np
from scipy.interpolate import interp1d, interp2d
import pickle
import os



rcp_efits = { # XRCP
             "XP184154_1":"184154_2940", "XP184154_2":"184154_4430",
             "XP184155_1":"184155_3000", # No second plunge.
             "XP184170_1":"184170_3000", "XP184170_2":"184170_4500",
             "XP184171_1":"184171_3000", "XP184171_2":"184171_4500",
             "XP184172_1":"184172_2500", "XP184172_2":"184172_4000",
             "XP184173_1":"184173_2500", "XP184173_2":"184173_4000",
             "XP184174_1":"184174_1600", "XP184174_2":"184174_3500",
             "XP184175_1":"184175_1600", "XP184175_2":"184175_3500",
             "XP184176_1":"184176_1600", "XP184176_2":"184176_3500",
             "XP184177_1":"184177_1600", "XP184177_2":"184177_3500",
             "XP184178_1":"184178_1600", "XP184178_2":"184178_3500",
             "XP184179_1":"184179_1600", "XP184179_2":"184179_3500",
             "XP184180_1":"184180_1600", "XP184180_2":"184180_3500",
             "XP184181_1":"184181_1600", "XP184181_2":"184181_3500",
             "XP184182_1":"184182_1600", "XP184182_2":"184182_3500",
             "XP184183_1":"184183_1600", # No second plunge.
             "XP184184_1":"184184_1600", "XP184184_2":"184184_3000",
             "XP184264_1":"184264_1530", "XP184264_2":"184264_2940",
             "XP184266_1":"184266_1530", "XP184266_2":"184266_3340",
             "XP184267_1":"184267_1600", "XP184267_2":"184264_5000",
             "XP184268_1":"184268_1600", "XP184268_2":"184268_3400",
             "XP184270_1":"184270_1600", "XP184270_2":"184270_3400",  # EFIT01. EFIT02 apparently didn't run correctly.
             "XP184271_1":"184271_1600", "XP184271_2":"184271_3400",
             "XP184272_1":"184272_1600", "XP184272_2":"184272_3400",
             "XP184273_1":"184273_1600", "XP184273_2":"184273_3400",
             "XP184274_1":"184274_1600", "XP184274_2":"184274_3400",
             "XP184527_1":"184527_1600", "XP184527_2":"184527_3400",
             "XP184528_1":"184528_1580", "XP184528_2":"184528_3380",
             "XP184529_1":"184529_1600", "XP184529_2":"184529_3400",
             "XP184530_1":"184530_1600", "XP184530_2":"184530_3400",
             "XP184531_1":"184531_1530", "XP184531_2":"184531_3340",
             "XP184532_1":"184532_1600", "XP184532_2":"184532_3400",
             "XP184533_1":"184533_1530", "XP184533_2":"184533_3340",
             "XP184535_1":"184535_1540", "XP184535_2":"184535_3340",
             "XP184536_1":"184536_1600", "XP184536_2":"184536_3400",
             "XP184537_1":"184537_1540", "XP184537_2":"184537_3340",
             "XP184538_1":"184538_1600", "XP184538_2":"184538_3400",

              # MRCP
             "MP184154_1":"184154_3000", "MP184154_2":"184154_4500",
             "MP184155_1":"184155_3000", # No second plunge.
             "MP184170_1":"184170_3000", "MP184170_2":"184170_4500",
             "MP184171_1":"184171_3000", "MP184171_2":"184171_4500",
             "MP184172_1":"184172_2500", "MP184172_2":"184172_4000",
             "MP184173_1":"184173_2500", "MP184173_2":"184173_4000",
             "MP184174_1":"184174_1600", "MP184174_2":"184174_3500",
             "MP184175_1":"184175_1600", "MP184175_2":"184175_3500",
             "MP184176_1":"184176_1600", "MP184176_2":"184176_3500",
             "MP184177_1":"184177_1600", "MP184177_2":"184177_3500",
             "MP184178_1":"184178_1600", "MP184178_2":"184178_3500",
             "MP184179_1":"184179_1600", "MP184179_2":"184179_3500",
             "MP184180_1":"184180_1600", # No second plunge.
             "MP184182_1":"184182_1620", "MP184182_2":"184182_3500",
             "MP184264_1":"184264_1600", "MP184264_2":"184264_3000",
             "MP184266_1":"184266_1600", "MP184266_2":"184266_3400",
             "MP184267_1":"184267_1610", "MP184267_2":"184267_5000",
             "MP184268_1":"184268_1600", "MP184268_2":"184268_3400",
             "MP184270_1":"184270_1600", "MP184270_2":"184270_3400",
             "MP184527_1":"184527_1600", "MP184527_2":"184527_3400",
             "MP184528_1":"184528_1600", "MP184528_2":"184528_3400",
             "MP184529_1":"184529_1600", "MP184529_2":"184529_3400",
             "MP184530_1":"184530_1610", "MP184530_2":"184530_3400",
             "MP184531_1":"184531_1600", "MP184531_2":"184531_3400",
             "MP184532_1":"184532_1610", "MP184532_2":"184532_3400",
             "MP184533_1":"184533_1600", "MP184533_2":"184533_3400",
             "MP187103_1":"187103_1620", "MP187103_2":"187103_3410",
             "MP187104_1":"187104_1620", "MP187104_2":"187104_3410",
             "MP187105_1":"187105_1620", "MP187105_2":"187105_3410",
             "MP187106_1":"187106_1620", "MP187106_2":"187106_3410",
             "MP187107_1":"187107_1620", "MP187107_2":"187107_3410",
             "MP187108_1":"187108_1620", "MP187108_2":"187108_3410",
             "MP187109_1":"187109_1620", "MP187109_2":"187109_3410",
             "MP187110_1":"187110_1620", "MP187110_2":"187110_3410",
             "MP187111_1":"187111_1620", "MP187111_2":"187111_3410"
             }

output_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"

# Load in every plunge as a DataFrame.
rcp_dfs = []; rcp_names = []
#root = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/"
root = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/"
filenames = os.listdir(root)
mp_names = ["Point", "Time (ms)", "R (cm)", "Rho", "Isat (A)", "ne (1e18 m-3)",
  "Te (eV)", "q_par (MW/m2)", "Mach", "Vflow (m/s)", "Vf1 (V)", "Vf2 (V)",
  "Vf3 (V)"]
xp_names = ["Point", "Time (ms)", "Z (cm)", "Rho", "Isat (A)", "ne (1e18 m-3)",
  "Te (eV)", "q_par (MW/m2)", "Mach", "Vflow (m/s)", "Vf (V)"]

for filename in np.sort(filenames):

    # Load differently so we can assign the correct column names.
    if filename[:2] == "MP":
        rcp_df = pd.read_csv(root + filename, sep="\t", names=mp_names, header=0)

    elif filename[:2] == "XP":
        rcp_df = pd.read_csv(root + filename, sep="\t", names=xp_names, header=0)

        # Need to fix a mistake by Dmitry here. He's given us distance from the
        # floor it seems. So let's add on the floor value to make it a true Z (m).
        floor_z = -1.244 * 100  # m to cm
        rcp_df["Z (cm)"] = rcp_df["Z (cm)"] + floor_z

    else:
        continue

    rcp_name = filename.split(".tab")[0]

    try:
        gfile_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/gfile_data_for_rcp/" + \
          "{}".format(rcp_efits[rcp_name])
    except KeyError:
        continue

    with open(gfile_path, "rb") as f:
        gfile = pickle.load(f)

    rcp_type = filename[:2]
    shot = int(filename[2:8])
    rcp_names.append(rcp_name)

    # Poloidal and toroidal fields.
    f_bp = interp2d(gfile["R"], gfile["Z"], gfile["Bp"])
    f_bt = interp2d(gfile["R"], gfile["Z"], gfile["Bt"])

    if rcp_type == "MP":
        r = rcp_df["R (cm)"] / 100
        z = np.full(len(rcp_df), -0.188)
    elif rcp_type == "XP":
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
    rcp_df["B"] = b
    rcp_df["Bt"] = bt
    rcp_df["Bp"] = bp

    # Calculate distance from separatrix. Limit to just portions of the LCFS
    # according to the RCP since it isn't a well-defined problem (fails the straight
    # line test).
    if rcp_type == "MP":
        mask = gfile["RBBBS"] > 2.00
        f_lcfs_r = interp1d(gfile["ZBBBS"][mask], gfile["RBBBS"][mask])
        rsep = f_lcfs_r(z[0])
        rcp_dist = r - rsep
    elif rcp_type == "XP":
        mask = gfile["ZBBBS"] < -0.80
        f_lcfs_z = interp1d(gfile["RBBBS"][mask], gfile["ZBBBS"][mask])
        zsep = f_lcfs_z(r[0])
        rcp_dist = zsep - z  # Could do abs, but do this so it's positive when we want it to be.

    f_psin = interp2d(gfile["R"], gfile["Z"], gfile["Psin"])

    if rcp_type == "MP":
        rcp_psin = f_psin(r, z)[0][::-1]
    elif rcp_type == "XP":
        rcp_psin = f_psin(r, z)[:, 0]

    # Additional mapping to map it all to R-Rsep OMP.
    # Need the R, Z of the magnetic axis which would require doing this all over
    # again lol.


    # Put it all into the rcp_df to make a new Excel file out of.
    rcp_df["Psin"] = rcp_psin
    rcp_df["Dist. from sep."] = rcp_dist
    rcp_dfs.append(rcp_df)

# Save to an Excel file with multiple sheets.
sort_idx  = np.argsort(rcp_names)
rcp_names = np.array(rcp_names)[sort_idx]
rcp_dfs   = np.array(rcp_dfs, dtype="object")[sort_idx]
with pd.ExcelWriter(output_path) as writer:
    for i in range(0, len(rcp_dfs)):
        rcp_dfs[i].to_excel(writer, sheet_name=rcp_names[i])
