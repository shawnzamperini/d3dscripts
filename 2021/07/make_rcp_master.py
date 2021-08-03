import pandas as pd
import os


output_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/rcp_master.xlsx"

# Load in every plunge as a DataFrame.
rcp_dfs = []; rcp_names = []
root = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/"
filenames = os.listdir(root)
names = ["Point", "Time (ms)", "R (cm)", "Rho", "Isat (A)", "ne (1e18 m-3)",
  "Te (eV)", "q_par (MW/m2)", "Mach", "Vflow (m/s)", "Vf1 (V)", "Vf2 (V)",
  "Vf3 (V)"]
for filename in filenames:
    if filename[:2] in ["MP", "XP"]:
        rcp_names.append(filename[:-4])
        rcp_dfs.append(pd.read_csv(root + filename, sep="\t", names=names, header=0))

# Save to an Excel file with multiple sheets.
with pd.ExcelWriter(output_path) as writer:
    for i in range(0, len(rcp_dfs)):
        rcp_dfs[i].to_excel(writer, sheet_name=rcp_names[i])
