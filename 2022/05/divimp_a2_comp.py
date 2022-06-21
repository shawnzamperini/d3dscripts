# Compare the crude collector probe model in DIVIMP to RBS data.
import oedge_plots
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# Get CP DIVIMP data.
def divimp_cp(op):
    cp = op.collector_probe(2.27, -0.188, 2.37, -0.188, showplot=False, numlocs=30)
    cp_itf_x = cp["r"]
    cp_itf_y = np.array(cp["flux1"])
    cp_otf_x = cp["r"]
    cp_otf_y = np.array(cp["flux2"])
    return {"itfx":cp_itf_x, "itfy":cp_itf_y, "otfx":cp_otf_x, "otfy":cp_otf_y}

root = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-div-"
cps = {}
for i in [1, 2, 3, 4, 5]:
    print("vscan{}".format(i))
    op = oedge_plots.OedgePlots(root+"vscan{}.nc".format(i))
    cps["vscan{}".format(i)] = divimp_cp(op)
for i in [1, 3, 4, 5, 6]:
    print("fscan{}".format(i))
    op = oedge_plots.OedgePlots(root+"fscan{}.nc".format(i))
    cps["fscan{}".format(i)] = divimp_cp(op)

# Get A2 data from sheet.
a2path = "/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
a2 = pd.read_excel(a2path)
a2_shift = 0.01
a2_itf_x = a2["R D (cm)"].values / 100 + a2_shift
a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values * 1e15 * 1e4
a2_otf_x = a2["R U (cm)"].values / 100 + a2_shift
a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values * 1e15 * 1e4

# Options
plot_type = "single"
exp_time = 4
absfac_mod = 0.2

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True, figsize=(8, 5))
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
if plot_type == "vscan":
    count = 0
    labels = [200, 400, 600, 800, 1000]
    for i in [1, 2, 3, 4, 5]:
        run = plot_type + str(i)
        cp_itf_x = cps[run]["itfx"]
        cp_itf_y = cps[run]["itfy"] * exp_time * absfac_mod
        cp_otf_x = cps[run]["otfx"]
        cp_otf_y = cps[run]["otfy"] * exp_time * absfac_mod
        ax1.plot(cp_itf_x, cp_itf_y, color=colors[count], label=labels[count])
        ax2.plot(cp_otf_x, cp_otf_y, color=colors[count], label=labels[count])
        count += 1
elif plot_type == "fscan":
    count = 0
    labels = [500, 1000, 1500, 2000, 2500, 3000]
    for i in [1, 3, 4, 5, 6]:
        run = plot_type + str(i)
        cp_itf_x = cps[run]["itfx"]
        cp_itf_y = cps[run]["itfy"] * exp_time * absfac_mod
        cp_otf_x = cps[run]["otfx"]
        cp_otf_y = cps[run]["otfy"] * exp_time * absfac_mod
        ax1.plot(cp_itf_x, cp_itf_y, color=colors[count], label=labels[count])
        ax2.plot(cp_otf_x, cp_otf_y, color=colors[count], label=labels[count])
        count += 1
elif plot_type == "single":
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-div11.nc"
    #ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-diffusion2.nc"
    op = oedge_plots.OedgePlots(ncpath)
    cp = op.collector_probe(2.27, -0.188, 2.37, -0.188, showplot=False, numlocs=30)
    cp_itf_x = cp["r"]
    cp_itf_y = np.array(cp["flux1"]) * exp_time * absfac_mod
    cp_otf_x = cp["r"]
    cp_otf_y = np.array(cp["flux2"]) * exp_time * absfac_mod
    ax1.plot(cp_itf_x, cp_itf_y, color="tab:purple")
    ax2.plot(cp_otf_x, cp_otf_y, color="tab:red")

ax1.scatter(a2_itf_x, a2_itf_y, color="tab:purple")
ax2.scatter(a2_otf_x, a2_otf_y, color="tab:red")
fig.supxlabel("R (m)")
ax1.set_ylabel("W Areal Density (m-2)")
ax1.set_title("ITF")
ax2.set_title("OTF")
ax1.set_xlim([2.28, 2.35])
ax1.set_ylim([0, 2e19])
ax2.legend()
fig.tight_layout()
fig.show()
