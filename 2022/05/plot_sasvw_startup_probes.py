import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


scan_rate = 500  # um/s
root = "/Users/zamperini/My Drive/Research/Data/lams_data/"
probes = ["SMCP06L", "SMCP06R", "SMCP07L", "SMCP07R", "SMCP08L", "SMCP08R", "DCP07L", "DCP07R"]
data = {}
for p in probes:
    df = pd.read_csv("{}{}.csv".format(root, p), skiprows=4, names=["time(s)", "w180", "w182", "w183", "w184", "w186"]).dropna()
    x = df["time(s)"].astype(float).values * (scan_rate / 10000)  # um to cm
    y = df[["w180", "w182", "w183", "w184", "w186"]].sum(axis=1).values

    # Convert LAMS counts to 1e17 atoms/cm2.
    y *= 2.23e6

    x = x - 1.6  # Offset

    ys = savgol_filter(y, 21, 2)

    data[p] = {"x":x, "y":y, "ys":ys}

# SMCP06L = OTF, SMCP07L, 08L = ITF
# DCP07L = OTF (I think)

# Will do a row of 3 plots.
lw = 2
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(11, 4), sharex=True, sharey=True)

ax1.set_title("SAS-VW Startup")
ax1.plot(data["SMCP06L"]["x"], data["SMCP06L"]["y"], label="OTF", color="tab:red", lw=lw)
ax1.plot(data["SMCP06R"]["x"], data["SMCP06R"]["y"], label="ITF", color="tab:purple", lw=lw)
ax1.grid()
ax1.legend(fontsize=12)

ax2.set_title("Background W (Upper Divertor SP)")
ax2.plot(data["SMCP07R"]["x"], data["SMCP07R"]["y"], label="OTF", color="tab:red", lw=lw)
ax2.plot(data["SMCP07L"]["x"], data["SMCP07L"]["y"], label="ITF", color="tab:purple", lw=lw)
ax2.grid()

ax3.set_title("DiMES W Pellet (Upper Divertor SP)")
ax3.plot(data["SMCP08R"]["x"], data["SMCP08R"]["y"], label="OTF", color="tab:red", lw=lw)
ax3.plot(data["SMCP08L"]["x"], data["SMCP08L"]["y"], label="ITF", color="tab:purple", lw=lw)
ax3.grid()

ax1.set_ylim([0, 5e10])
ax1.set_xlim([0, 8])
fig.supxlabel("Distance along probe (cm)", fontsize=14)
ax1.set_ylabel(r"W Areal Density W/$\mathdefault{cm^2}$", fontsize=14)

fig.tight_layout()
fig.show()

fig, ax1 = plt.subplots()
ax1.set_title("SAS-VW Startup (DiMES)")
ax1.plot(data["DCP07R"]["x"], data["DCP07R"]["y"], label="OTF", color="tab:red", lw=lw)
ax1.plot(data["DCP07L"]["x"], data["DCP07L"]["y"], label="ITF", color="tab:purple", lw=lw)
ax1.grid()
ax1.legend(fontsize=12)
ax1.set_xlabel("Distance along probe (cm)", fontsize=14)
ax1.set_ylabel(r"W Areal Density W/$\mathdefault{cm^2}$", fontsize=14)
ax1.set_ylim([0, 2e12])
ax1.set_xlim([0, 4])
fig.tight_layout()
fig.show()
