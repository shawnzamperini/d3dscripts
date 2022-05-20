import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

root = "/Users/zamperini/My Drive/Research/Data/lams_data/"

# For all these, L = ITF and R = OTF
labels = [187113, 187116, 187119, 187122]
cp1l = pd.read_excel(root+"SMCP01L2_Analysis.xlsx")  # 187113
cp1r = pd.read_excel(root+"SMCP01R2_Analysis.xlsx")
cp2l = pd.read_excel(root+"SMCP02L2_Analysis.xlsx")  # 187116
cp2r = pd.read_excel(root+"SMCP02R2_Analysis.xlsx")
cp3l = pd.read_excel(root+"SMCP03L2_Analysis.xlsx")  # 187119
cp3r = pd.read_excel(root+"SMCP03R2_Analysis.xlsx")
cp4l = pd.read_excel(root+"SMCP04L2_Analysis.xlsx")  # 187122
cp4r = pd.read_excel(root+"SMCP04R2_Analysis.xlsx")


def get_xy_c13(df_itf, df_otf):

    itfx = df_itf["Axial Location [mm]"] / 10
    itfy = df_itf["13C Excess"]
    itfy_si = df_itf["Total Si"]
    otfx = df_otf["Axial Location [mm]"] / 10
    otfy = df_otf["13C Excess"]
    otfy_si = df_otf["Total Si"]

    # Convert LAMS counts to 1e17 atoms/cm2.
    itfy = (itfy + 346.05) / 11942
    otfy = (otfy + 346.05) / 11942
    itfy = itfy * 1e17
    otfy = otfy * 1e17

    return itfx, itfy, itfy_si, otfx, otfy, otfy_si

itfx1, itfy1, itfy1_si, otfx1, otfy1, otfy1_si = get_xy_c13(cp1l, cp1r)
itfx2, itfy2, itfy2_si, otfx2, otfy2, otfy2_si = get_xy_c13(cp2l, cp2r)
itfx3, itfy3, itfy3_si, otfx3, otfy3, otfy3_si = get_xy_c13(cp3l, cp3r)
itfx4, itfy4, itfy4_si, otfx4, otfy4, otfy4_si = get_xy_c13(cp4l, cp4r)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 6))

ax1.plot(itfx1, itfy1, label=187113)
ax1.plot(itfx2, itfy2, label=187116)
ax1.plot(itfx3, itfy3, label=187119)
ax1.plot(itfx4, itfy4, label=187122)
ax1.set_ylim([0, 1e17])

ax2.plot(otfx1, otfy1, label=187113)
ax2.plot(otfx2, otfy2, label=187116)
ax2.plot(otfx3, otfy3, label=187119)
ax2.plot(otfx4, otfy4, label=187122)
ax2.set_ylim([0, 1e17])
ax2.legend()

ax3.plot(itfx1, itfy1_si, label=187113)
ax3.plot(itfx2, itfy2_si, label=187116)
ax3.plot(itfx3, itfy3_si, label=187119)
ax3.plot(itfx4, itfy4_si, label=187122)
ax3.set_ylim([0, 2e7])

ax4.plot(otfx1, otfy1_si, label=187113)
ax4.plot(otfx2, otfy2_si, label=187116)
ax4.plot(otfx3, otfy3_si, label=187119)
ax4.plot(otfx4, otfy4_si, label=187122)
ax4.set_ylim([0, 2e7])

fig.tight_layout()
fig.show()
