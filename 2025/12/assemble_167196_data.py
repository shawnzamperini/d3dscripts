import numpy as np
import pandas as pd
import pickle
from ThomsonClass import ThomsonClass


exp_data = {}

# Load TS data
#ts_path = "/home/zamp/data/ts_167195.pickle"
#with open(ts_path, "rb") as f:
#	ts = pickle.load(f)

# Load data from core TS system
ts_core = ThomsonClass(167196, "core")
ts_core.load_ts(tunnel=False)
r = ts_core.ts_dict["r"]["Y"]
z = ts_core.ts_dict["z"]["Y"]
t = ts_core.ts_dict["time"]["Y"]

# Limit between a reasonable time range
tmin = 2500; tmax = 5000
time_mask = np.logical_and(t >= tmin, t <= tmax)
ts_core_dict = {}
ts_core_dict["shot"] = 167196
ts_core_dict["t_min (ms)"] = tmin
ts_core_dict["t_max (ms)"] = tmax
ts_core_dict["r (m)"] = r
ts_core_dict["z (m)"] = z
ts_core_dict["t (ms)"] = t[time_mask]
ts_core_dict["ne (m-3)"] = ts_core.ts_dict["density"]["Y"][:, time_mask]
ts_core_dict["te (eV)"] = ts_core.ts_dict["temp"]["Y"][:, time_mask]

# Put into a single dict
exp_data["ts_core"] = ts_core_dict

# Load data from divertor TS system. Use 167195 because strike point was
# sweeping and it fills out a 2D profile
#ts_core = ThomsonClass(167195, "divertor")
#ts_core.load_ts(tunnel=False)

# Load RCP data (#167192-5)
rcps = {}
cas = {}
rcp_Z = -0.188
for shot in range(167192, 167196):
	for plunge in [1, 2]:

		# Load into a DataFrame
		rcp_path = "/home/zamp/data/MP{}_{}.tab".format(shot, plunge)
		rcp = pd.read_csv(rcp_path, delimiter='\t')
		rcp["Z(m)"] = np.full(len(rcp), rcp_Z)

		# Better labels
		# ...

		# Plasma potential
		rcp["Vp (V)"] = (rcp["Vf1(V)"] + rcp["Vf2(V)"]) / 2 + 3 * rcp["Te(eV)"]
		
		# Convert to a dictionary (list tells it to not include the index)
		rcp_dict = rcp.to_dict("list")
		rcps["MP{}_{}".format(shot, plunge)] = rcp_dict

		# Similar treatment but for the conditional averaged data
		rcp_path = "/home/zamp/data/CA_{}_{}.tab".format(shot, plunge)
		try:
			rcp = pd.read_csv(rcp_path, delimiter='\t')

		# Only CA data for 167193 and 167195
		except:
			continue

		# Similar to the above
		rcp["Z(m)"] = np.full(len(rcp), rcp_Z)
		rcp_dict = rcp.to_dict("list")
		cas["CA{}_{}".format(shot, plunge)] = rcp_dict

# Put into dictionary
exp_data["rcp"] = rcps
exp_data["rcp_ca"] = cas

# Save as pickled dictionary
with open("/home/zamp/data/exp_data_167196.pickle", "wb") as f:
	pickle.dump(exp_data, f)
