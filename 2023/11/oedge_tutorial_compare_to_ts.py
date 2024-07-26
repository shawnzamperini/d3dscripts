import oedge_plots
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load Thomson scattering data.
corets_shift = 0.0
ts_path = "/Users/zamperini/Documents/d3d_work/divimp_files/oedge_tutorial/ts_167195.pickle"
with open(ts_path, "rb") as f:
    ts = pickle.load(f)
ts_plot = {"core": {}, "divertor": {}, "tangential": {}}
for sys in ts.keys():
    tmp = ts[sys]
    mask = np.logical_and(tmp["time"] >= 2500, tmp["time"] <= 5000)
    ts_plot[sys]["time"] = tmp["time"][mask]
    for key in ["te", "te_err", "ne", "ne_err", "psin"]:
        if sys == "core" and key == "psin":
            ts_plot[sys][key] = tmp[key][:, mask] + corets_shift
        else:
            ts_plot[sys][key] = tmp[key][:, mask]
    ts_plot[sys]["chord"] = tmp["chord"]

# Load the RCP data. Data has already been shifted inward by 1.5 cm due to EFIT uncertainties.
rcp_path = "/Users/zamperini/Documents/d3d_work/divimp_files/oedge_tutorial/rcp_156195_2.csv"
rcp = pd.read_csv(rcp_path)

# Load OEDGE run and extract a series of profiles along the locations of TS and RCP.
op_path = "/Users/zamperini/Documents/d3d_work/divimp_files/oedge_tutorial/d3d-167196-osm-v1.nc"
op = oedge_plots.OedgePlots(op_path)
op_tsc_te = op.along_line(1.94, 1.94, 0.67, 0.85, "KTEBS", "psin")
op_tsc_ne = op.along_line(1.94, 1.94, 0.67, 0.85, "KNBS", "psin")
op_tsd_te = op.along_line(1.484, 1.484, -0.82, -1.17, "KTEBS", "psin")
op_tsd_ne = op.along_line(1.484, 1.484, -0.82, -1.17, "KNBS", "psin")
op_rcp_te = op.along_line(2.18, 2.30, -0.188, -0.188, "KTEBS", "psin")
op_rcp_ne = op.along_line(2.18, 2.30, -0.188, -0.188, "KNBS", "psin")

# Switch to include ring numbers above every data point. Messy, but useful for zooming in and debugging.
labels = True

# Now we do our comparison plots.
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(8, 5))

# Core TS Te.
x = ts_plot["core"]["psin"].flatten()
y = ts_plot["core"]["te"].flatten()
yerr = ts_plot["core"]["te_err"].flatten()
ax1.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax1.plot(op_tsc_te["psin"], op_tsc_te["KTEBS"], color="tab:red")
ax1.set_xlabel("Psin")
ax1.set_title("Core TS Te")
ax1.set_xlim([0.99, 1.15])
ax1.set_ylim([0, 100])

# Core TS ne.
x = ts_plot["core"]["psin"].flatten()
y = ts_plot["core"]["ne"].flatten()
yerr = ts_plot["core"]["ne_err"].flatten()
ax4.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax4.plot(op_tsc_ne["psin"], op_tsc_ne["KNBS"], color="tab:red")
ax4.set_xlabel("Psin")
ax4.set_title("Core TS ne")
ax4.set_xlim([0.99, 1.15])
ax4.set_ylim([0, 2.0e19])

# Divertor TS Te
x = ts_plot["divertor"]["psin"].flatten()
y = ts_plot["divertor"]["te"].flatten()
yerr = ts_plot["divertor"]["te_err"].flatten()
ax2.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax2.plot(op_tsd_te["psin"], op_tsd_te["KTEBS"], color="tab:red")
ax2.set_xlabel("Psin")
ax2.set_title("Divertor TS Te")
ax2.set_xlim([0.99, 1.03])
ax2.set_ylim([0, 100])

# Divertor TS ne
x = ts_plot["divertor"]["psin"].flatten()
y = ts_plot["divertor"]["ne"].flatten()
yerr = ts_plot["divertor"]["ne_err"].flatten()
ax5.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax5.plot(op_tsd_ne["psin"], op_tsd_ne["KNBS"], color="tab:red")
ax5.set_xlabel("Psin")
ax5.set_title("Divertor TS ne")
ax5.set_xlim([0.99, 1.03])
ax5.set_ylim([0, 1e20])

# RCP Te.
x = rcp["psin"].values
y = rcp["Te(eV)"].values
ax3.scatter(x, y, s=15, color="k")
ax3.plot(op_rcp_te["psin"], op_rcp_te["KTEBS"], color="tab:red", marker=".")
ax3.set_xlabel("Psin")
ax3.set_title("RCP Te")
# ax3.axvline(2.2367, color="k", linestyle="--")
ax3.set_xlim([0.99, 1.3])
ax3.set_ylim([0, 50])

# RCP ne.
x = rcp["psin"].values
y = rcp["Ne(E18 m-3)"].values * 1e18
ax6.scatter(x, y, s=15, color="k")
ax6.plot(op_rcp_ne["psin"], op_rcp_ne["KNBS"], color="tab:red", marker=".")
ax6.set_xlabel("Psin")
ax6.set_title("RCP ne")
# ax6.axvline(2.2367, color="k", linestyle="--")
ax6.set_xlim([0.99, 1.3])
ax6.set_ylim([0, 2e19])

fig.tight_layout()
fig.show()

"""
  .. code-block:: python

    import pickle
    from os.path import expanduser
    
    # Set to True if you want to filter out ELMs.
    filter_elms = False
    
    # Load lists containing the raw data points.
    core = OMFIT["OMFITprofiles"]["OUTPUTS"]["RAW"]["TS"]["core_r+1"]
    div = OMFIT["OMFITprofiles"]["OUTPUTS"]["RAW"]["TS"]["divertor_r-1"]
    
    # Convert data into dictionaries with only the data we need. Core first.
    raw_data = {}
    for sys in [core, div]:
        for channel in range(0, len(sys)):
            c = sys[channel]
        
            # If desired, create mask filtering ELMs.
            if filter_elms:
                mask = [True if i==1 else False for i in c["selected_times"]]
            else:
                mask = [True for i in c["selected_times"]]
        
            time = c["time"].data[mask]
            r = c["R"].data[0][mask]
            z = c["Z"].data[0][mask]
            psin = c["psi_n"][0][mask]
            ne = c["n_e"][0][mask]
            te  = c["T_e"][0][mask]
    
            if sys == core:
                label = "{}_{}".format("core", channel)
            elif sys == div:
                label = "{}_{}".format("divertor", channel)
            
            raw_data[label] = {"time":time, "R":r, "Z":z, "psin":psin, "ne":ne, "Te":te}
    
    home = expanduser("~")
    fname = "{}/ts_dict.pickle".format(home)
    with open(fname, "wb") as f:
        pickle.dump(raw_data, f)
"""
