from gadata import gadata
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d


unfavorable = False
if unfavorable:
    diag = 184527
    coll = 184535
    dens_lims = [1.5e19, 3.5e19]
    rvs_lims = [1.30, 1.35]
    a_coords = (0.03, 0.55)
else:
    diag = 184267
    coll = 184271
    dens_lims = [1.5e19, 2.5e19]
    rvs_lims = [1.32, 1.36]
    a_coords = (0.03, 0.85)

tags = ["DENSV2", "PINJ", "UOB", "ASDEX_UPBAF", "FS04", "RVSOUT"]
diags = {}
colls = {}
for tag in tags:
    gadiag = gadata(tag, diag)
    gacoll = gadata(tag, coll)
    diags[tag] = [gadiag.xdata/1000, gadiag.zdata]
    colls[tag] = [gacoll.xdata/1000, gacoll.zdata]

# Fix weird shit with IP.
#signals247["IP"] = [signals247["IP"][0].value, signals247["IP"][1].value]
#signals277["IP"] = [signals277["IP"][0].value, signals277["IP"][1].value]
diags["UOB"] = [diags["UOB"][0].value, diags["UOB"][1].value]
colls["UOB"] = [colls["UOB"][0].value, colls["UOB"][1].value]

# Smooth PINJ.
diags["PINJ"][1] = savgol_filter(diags["PINJ"][1], 3001, 2)
colls["PINJ"][1] = savgol_filter(colls["PINJ"][1], 3001, 2)
diags["RVSOUT"][1] = savgol_filter(diags["RVSOUT"][1], 31, 2)
colls["RVSOUT"][1] = savgol_filter(colls["RVSOUT"][1], 31, 2)

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(8, 5), sharex=True)

ax1.plot(diags["DENSV2"][0], diags["DENSV2"][1], color="tab:red", label=diag)
ax1.plot(colls["DENSV2"][0], colls["DENSV2"][1], color="tab:purple", label=coll)
ax1.set_ylabel("Density (m-3)")
ax1.set_xlim([0, 6.5])
ax1.set_ylim(dens_lims)
ax1.legend(ncol=2)
ax1.text(a_coords[0], a_coords[1], "a)", transform=ax1.transAxes)

ax3.plot(diags["PINJ"][0], diags["PINJ"][1]/1000, color="tab:red")
ax3.plot(colls["PINJ"][0], colls["PINJ"][1]/1000, color="tab:purple")
ax3.set_ylabel("Pinj (MW)")
ax3.text(0.03, 0.85, "c)", transform=ax3.transAxes)

ax2.plot(diags["UOB"][0], diags["UOB"][1], color="tab:red")
ax2.plot(colls["UOB"][0], colls["UOB"][1], color="tab:purple")
ax2.set_ylabel("UOB (V)")
ax2.text(0.03, 0.85, "b)", transform=ax2.transAxes)

ax4.plot(diags["ASDEX_UPBAF"][0], diags["ASDEX_UPBAF"][1]*1000, color="tab:red")
ax4.plot(colls["ASDEX_UPBAF"][0], colls["ASDEX_UPBAF"][1]*1000, color="tab:purple")
ax4.set_ylabel("UOB Pressure (mTorr)")
ax4.text(0.03, 0.85, "d)", transform=ax4.transAxes)

ax6.plot(diags["FS04"][0], diags["FS04"][1], color="tab:red")
ax6.plot(colls["FS04"][0], colls["FS04"][1], color="tab:purple")
ax6.set_ylabel("Filterscope\n(ph/sr " + r"cm$^2$" + " s)")
ax6.set_ylim([0, 2e15])
ax6.set_xlabel("Time (s)")
ax6.text(0.03, 0.85, "f)", transform=ax6.transAxes)

ax5.plot(diags["RVSOUT"][0], diags["RVSOUT"][1], color="tab:red")
ax5.plot(colls["RVSOUT"][0], colls["RVSOUT"][1], color="tab:purple")
ax5.set_ylabel("Strike Point (m)")
ax5.set_ylim(rvs_lims)
ax5.set_xlabel("Time (s)")
ax5.text(0.03, 0.85, "e)", transform=ax5.transAxes)

fig.tight_layout()
fig.show()
