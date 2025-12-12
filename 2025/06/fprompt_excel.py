from package_redeposition_model.semi_analytical_redeposition_model import RedepositionModel
import pandas as pd
import numpy as np
from tqdm import tqdm


# Inputs
xlpath = "/mnt/g/My Drive/Research/Documents/2025/excel_wall_erosion_v2.xlsx"
target_mat = "W"  # W or Fe
sput_imp = "Ne"   # Ne, ...
#sput_charge = 1

print("Calculating results for {} on {}".format(sput_imp, target_mat))

# Load excel sheet
#  sudo mount -t drvfs G: /mnt/g
if target_mat == "W":
	sheet_name = "fnprompt_w"
elif target_mat == "Fe":
	sheet_name = "fnprompt_fe"
df = pd.read_excel(xlpath, sheet_name=sheet_name)

# Max charge for each element
max_charge = {"Ne":10}

# Loop through each combo of inputs
for charge in range(1, max_charge[sput_imp] + 1):
	col_name = sput_imp + str(charge)
	print(col_name)

	# Initialize model
	model = RedepositionModel(target_species=target_mat, imp_species=sput_imp, 
		qp=charge)

	# Calculate non-prompt redeposition fractions
	fnps = np.zeros(len(df))
	for i in tqdm(range(len(df))):
		results = model.run_model(TeSE=df["te"].iloc[i], neSE=df["ne"].iloc[i], 
			B=df["b"].iloc[i], alphaB=np.radians(df["ang (deg)"].iloc[i]))
		fnps[i] = results["fnp"]
	
	# Put into df
	df[col_name] = fnps

# Save df as a csv so one can just load it up and manually copy things over
fname = "{}_on_{}.csv".format(sput_imp, target_mat)
df.to_csv(fname)
print("Saved as: {}".format(fname))

# plasma background
#TeSE=35                  # eV
#neSE=1.2e19              # m^-3
#B=2                      # T
#alphaB= 3.14 / 180 * 1   # rad

# Run model
#results = model.run_model(TeSE = TeSE, neSE = neSE, B = B, alphaB = alphaB)
#print(results)
