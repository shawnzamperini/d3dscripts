import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


lams_path = "/Users/zamperini/My Drive/Research/Data/lams_data/methane_lams_master.xlsx"
df = pd.read_excel(lams_path, sheet_name="For Export")
raw_unf_itf_x = df["ml04_loc"].values
raw_unf_itf_c13 = df["ml04_excess_c13"].values
raw_unf_itf_si = df["ml04_tot_si"].values
raw_unf_otf_x = df["mr04_loc"].values
raw_unf_otf_c13 = df["mr04_excess_c13"].values
raw_unf_otf_si = df["mr04_tot_si"].values

raw_fav_itf_x = df["mr21_loc"].values
raw_fav_itf_c13 = df["mr21_excess_c13"].values
raw_fav_itf_si = df["mr21_tot_si"].values
raw_fav_otf_x = df["ml21_loc"].values
raw_fav_otf_c13 = df["ml21_excess_c13"].values
raw_fav_otf_si = df["ml21_tot_si"].values

# Remove data before zero, smoothing, background subtraction.
def process(lams_itf_x, lams_itf_y, lams_otf_x, lams_otf_y):
    mask1 = lams_itf_x > 0
    mask2 = lams_otf_x > 0
    lams_itf_x = lams_itf_x[mask1]
    lams_itf_y = lams_itf_y[mask1]
    lams_otf_x = lams_otf_x[mask2]
    lams_otf_y = lams_otf_y[mask2]
    lams_itf_ys = savgol_filter(lams_itf_y, 51, 2)
    lams_otf_ys = savgol_filter(lams_otf_y, 51, 2)
    min_itfy = lams_itf_ys.min()
    min_otfy = lams_otf_ys.min()
    lams_itf_y = lams_itf_y - min_itfy
    lams_otf_y = lams_otf_y - min_otfy
    lams_itf_ys = lams_itf_ys - min_itfy
    lams_otf_ys = lams_otf_ys - min_otfy

    return lams_itf_x, lams_itf_ys, lams_otf_x, lams_otf_ys

unf_itf_x, unf_itf_c13, unf_otf_x, unf_otf_c13 = process(raw_unf_itf_x, raw_unf_itf_c13, raw_unf_otf_x, raw_unf_otf_c13)
unf_itf_x, unf_itf_si, unf_otf_x, unf_otf_si = process(raw_unf_itf_x, raw_unf_itf_si, raw_unf_otf_x, raw_unf_otf_si)
fav_itf_x, fav_itf_c13, fav_otf_x, fav_otf_c13 = process(raw_fav_itf_x, raw_fav_itf_c13, raw_fav_otf_x, raw_fav_otf_c13)
fav_itf_x, fav_itf_si, fav_otf_x, fav_otf_si = process(raw_fav_itf_x, raw_fav_itf_si, raw_fav_otf_x, raw_fav_otf_si)

fontsize = 14

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 5), sharex=True)

ax1.plot(unf_itf_x, unf_itf_c13, color="tab:red")
ax1.set_ylim([0, 1750])
ax1.set_xlim([0, 10])
ax1.set_ylabel("ITF C13 Counts", color="tab:red", fontsize=fontsize)
ax1.tick_params(axis="y", color="tab:red", labelcolor="tab:red")
ax1.set_title("#184527", fontsize=fontsize)
ax11 = ax1.twinx()
ax11.plot(unf_itf_x, unf_itf_si, color="tab:purple")
ax11.set_ylim([0, 1.2e6])
ax11.tick_params(axis="y", color="tab:purple", labelcolor="tab:purple")

ax3.plot(unf_otf_x, unf_otf_c13, color="tab:red")
ax3.set_ylim([0, 1750])
ax3.tick_params(axis="y", color="tab:red", labelcolor="tab:red")
ax3.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
ax3.set_ylabel("OTF C13 Counts", color="tab:red", fontsize=fontsize)
ax33 = ax3.twinx()
ax33.plot(unf_otf_x, unf_otf_si, color="tab:purple")
ax33.set_ylim([0, 1.2e6])
ax33.tick_params(axis="y", color="tab:purple", labelcolor="tab:purple")

ax2.plot(fav_itf_x, fav_itf_c13, color="tab:red", label="C13")
ax2.set_ylim([0, 1750])
ax2.tick_params(axis="y", color="tab:red", labelcolor="tab:red")
ax2.set_title("#184267", fontsize=fontsize)
ax22 = ax2.twinx()
ax22.plot(fav_itf_x, fav_itf_si, color="tab:purple", label="Si")
ax22.set_ylim([0, 1.2e6])
ax22.set_ylabel("ITF Si Counts", color="tab:purple", fontsize=fontsize)
ax22.tick_params(axis="y", color="tab:purple", labelcolor="tab:purple")

ax4.plot(fav_otf_x, fav_otf_c13, color="tab:red")
ax4.set_ylim([0, 1750])
ax4.tick_params(axis="y", color="tab:red", labelcolor="tab:red")
ax4.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
ax44 = ax4.twinx()
ax44.plot(fav_otf_x, fav_otf_si, color="tab:purple")
ax44.set_ylim([0, 1.2e6])
ax44.tick_params(axis="y", color="tab:purple", labelcolor="tab:purple")
ax44.set_ylabel("OTF Si Counts", color="tab:purple", fontsize=fontsize)

fig.tight_layout()
fig.show()
