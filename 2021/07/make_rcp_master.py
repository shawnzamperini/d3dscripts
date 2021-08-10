import pandas as pd
import os


#output_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/rcp_master.xlsx"
output_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/rcp_master.xlsx"

# Load in every plunge as a DataFrame.
rcp_dfs = []; rcp_names = []
#root = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/"
root = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/"
filenames = os.listdir(root)
mp_names = ["Point", "Time (ms)", "R (cm)", "Rho", "Isat (A)", "ne (1e18 m-3)",
  "Te (eV)", "q_par (MW/m2)", "Mach", "Vflow (m/s)", "Vf1 (V)", "Vf2 (V)",
  "Vf3 (V)"]
xp_names = ["Point", "Time (ms)", "Z (cm)", "Rho", "Isat (A)", "ne (1e18 m-3)",
  "Te (eV)", "q_par (MW/m2)", "Mach", "Vflow (m/s)", "Vf"]
for filename in filenames:
    if filename[:2] == "MP":
        rcp_names.append(filename[:-4])
        rcp_dfs.append(pd.read_csv(root + filename, sep="\t", names=mp_names, header=0))
    if filename[:2] == "XP":
        rcp_names.append(filename[:-4])
        xp_df = pd.read_csv(root + filename, sep="\t", names=xp_names, header=0)

        # Need to fix a mistake by Dmitry here. He's given us distance from the
        # floor it seems. So let's add on the floor value to make it a true Z (m).
        floor_z = -1.244 * 100  # m to cm
        xp_df["Z (cm)"] = xp_df["Z (cm)"] + floor_z
        rcp_dfs.append(xp_df)


# Save to an Excel file with multiple sheets.
with pd.ExcelWriter(output_path) as writer:
    for i in range(0, len(rcp_dfs)):
        rcp_dfs[i].to_excel(writer, sheet_name=rcp_names[i])
