import MDSplus
import pandas as pd
import numpy  as np
from   gadata import gadata
import matplotlib.pyplot as plt


unf_shots = np.append(np.arange(184152, 184164), np.arange(184524, 184541))
fav_shots = np.append(np.arange(184164, 184185), np.append(np.arange(184262, 184275), np.arange(187103, 187112)))
all_shots = np.append(unf_shots, fav_shots)
conn = MDSplus.Connection("localhost")

rcp_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/rcp_master.xlsx"
dfs = pd.read_excel(rcp_path, sheet_name=None)

def avg_val(shot, tag, time, time_width=10):
    gaobj = gadata(tag, shot, connection=conn, print_out=False)
    keep = np.logical_and(gaobj.xdata>=time-time_width, gaobj.xdata<=time+time_width)
    avg = gaobj.zdata[keep].mean()
    return avg

meth_df = pd.DataFrame()
#for shot in all_shots:
for rcp_name in dfs.keys():

    shot = int(rcp_name[2:8])
    rcp = dfs[rcp_name]

    print(shot)

    t = rcp["Time(ms)"].mean()
    s = pd.Series(name=str(shot)+"_{:d}".format(int(t)), dtype=np.float)
    #s2 = pd.Series(name=str(shot)+"_2", dtype=np.float)

    mask1 = np.logical_and(rcp["Rho"] >= 1.00, rcp["Rho"] < 1.02)
    mask2 = np.logical_and(rcp["Rho"] >= 1.02, rcp["Rho"] < 1.04)
    mask3 = np.logical_and(rcp["Rho"] >= 1.04, rcp["Rho"] < 1.06)
    mask4 = rcp["Rho"] >= 1.06
    s["mach_10-12"] = rcp["Machn"][mask1].mean()
    s["mach_12-14"] = rcp["Machn"][mask2].mean()
    s["mach_14-16"] = rcp["Machn"][mask3].mean()
    s["mach_16-"] = rcp["Machn"][mask4].mean()

    # Load relevant info, averaging at each time a plunge was (about 1600 and 3500).
    #for s in [s1, s2]:
    if shot in unf_shots:
        s["direction"] = "unfavorable"
    else:
        s["direction"] = "favorable"
    s["ne"] = avg_val(shot, "DENSV2", t)
    s["pinj"] = avg_val(shot, "PINJ", t) / 1000  # kW to MW
    s["rvsout"] = avg_val(shot, "RVSOUT", t)
    s["ip"] = avg_val(shot, "IP", t)
    s["prad_core"] = avg_val(shot, "PRAD_CORE", t) / 1e6  # W to MW
    s["psol"] = s["pinj"] - s["prad_core"]
    s["wmhd"] = avg_val(shot, "WMHD", t)
    s["gapbot"] = avg_val(shot, "GAPBOT", t)
    s["aminor"] = avg_val(shot, "AMINOR", t)
    s["f_greenwald"] = s["ne"] / (s["ip"] * 10**(-6) / (np.pi * s["aminor"]**2) * 1e20)
    s["time"] = t
    s["pfile"] = rcp_name
    s["ptype"] = rcp_name[:2]

    # Integrate GASA to determine total injected gas.
    gaobj_gasa = gadata("GASA", shot, connection=conn, print_out=False)
    gasa_int = np.trapz(gaobj_gasa.zdata, gaobj_gasa.xdata)
    s["total_gasa"] = gasa_int

    # Put into master DataFrame.
    meth_df = meth_df.append(s)

meth_df.to_excel("methane_shots.xlsx")
