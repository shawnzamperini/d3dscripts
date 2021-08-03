import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# B6 requires manual handling.
plist = ["B4", "B5", "B7", "B9", "B10",
         "C4", "C5", "C6", "C7", "C8", "C10"]
plist_rev = ["B3", "B4", "B5", "C3", "C4", "C5"]
plist_for = ["B6", "B7", "B8", "B9", "B10", "C6", "C7", "C8", "C9", "C10"]

# Some probe don't have RBS, so to get Total ITF/OTF we have to use LAMS.
plist_only_lams = ["B4", "C4", "B6", "B9", "C6"]

rbs_root = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/"
lams_root = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/"

# Lists to hold the results.
rev_all_tot_itfotf      = []
rev_all_tot_itfotf_err  = []
rev_all_peak_otfitf     = []
rev_all_peak_otfitf_err = []
for_all_tot_itfotf      = []
for_all_tot_itfotf_err  = []
for_all_peak_otfitf     = []
for_all_peak_otfitf_err = []

# Data from the set of 3DLIM runs.
lim_otfitf_peak = [1.36, 1.07, 0.99, 0.94, 0.72]
lim_itfotf_tot  = [1.91, 1.08, 1.01, 0.94, 0.50]

for p in plist:

    print(p)

    # Load in RBS data if it's there.
    if p in plist_only_lams:
        pass
    else:
        rbs_df = pd.read_excel("{}{}.xlsx".format(rbs_root, p))

    # Load in LAMS data, padding with zeros if needed and assigning correct
    # side as ITF or OTF.
    if p[0] == "B":
        if int(p.split("B")[1]) < 10:
            p_paddedD = "BD0" + p.split("B")[1]
            p_paddedU = "BU0" + p.split("B")[1]
        else:
            p_paddedD = "BD" + p.split("B")[1]
            p_paddedU = "BU" + p.split("B")[1]
    else:
        if int(p.split("C")[1]) < 10:
            p_paddedD = "CD0" + p.split("C")[1]
            p_paddedU = "CU0" + p.split("C")[1]
        else:
            p_paddedD = "CD" + p.split("C")[1]
            p_paddedU = "CU" + p.split("C")[1]

    if p in plist_rev:
        itf_df = pd.read_excel("{}{}_Map_Analysis.xlsx".format(lams_root,
                 p_paddedD)).pivot_table(columns='z Location [mm]',
                 index='Axial Location [mm]')
        otf_df = pd.read_excel("{}{}_Map_Analysis.xlsx".format(lams_root,
                 p_paddedU)).pivot_table(columns='z Location [mm]',
                 index='Axial Location [mm]')
    else:
        itf_df = pd.read_excel("{}{}_Map_Analysis.xlsx".format(lams_root,
                 p_paddedU)).pivot_table(columns='z Location [mm]',
                 index='Axial Location [mm]')
        otf_df = pd.read_excel("{}{}_Map_Analysis.xlsx".format(lams_root,
                 p_paddedD)).pivot_table(columns='z Location [mm]',
                 index='Axial Location [mm]')

    # Calculate total W on ITF and OTF sides within first 5 cm.
    rad_cutoff = 5.0
    if p in plist_only_lams:
        pass
    else:
        sum_range = rbs_df["Distance from Tip D (cm)"] <= rad_cutoff
        if p in plist_rev:
            itf_side = "D"
            otf_side = "U"
        else:
            itf_side = "U"
            otf_side = "D"

        # Get relevant arrays and values from dataframes.
        den_itf = rbs_df["W Areal Density {} (1e15 W/cm2)".format(itf_side)][sum_range].values
        den_otf = rbs_df["W Areal Density {} (1e15 W/cm2)".format(otf_side)][sum_range].values
        tot_itf = rbs_df["W Areal Density {} (1e15 W/cm2)".format(itf_side)][sum_range].sum()
        tot_otf = rbs_df["W Areal Density {} (1e15 W/cm2)".format(otf_side)][sum_range].sum()
        err_itf = rbs_df["W Areal Density Error {} (1e15 W/cm2)".format(itf_side)][sum_range].values
        err_otf = rbs_df["W Areal Density Error {} (1e15 W/cm2)".format(otf_side)][sum_range].values

        # Calculate ratio and propagate errors.
        tot_itf_err = np.sqrt(np.square(err_itf).sum())
        tot_otf_err = np.sqrt(np.square(err_otf).sum())
        tot_itfotf = tot_itf / tot_otf
        tot_itfotf_err = tot_itfotf * np.sqrt(np.square(tot_itf_err / tot_itf) + np.square(tot_otf_err / tot_otf))
        print("  Total ITF/OTF: {:.2f} +/- {:.2f}".format(tot_itfotf, tot_itfotf_err))

    # Calculate the average (normalized) poloidal profiles.
    zloc_otf = otf_df.iloc[otf_df.index<=rad_cutoff*10].mean()['Total W'].index
    zloc_itf = itf_df.iloc[itf_df.index<=rad_cutoff*10].mean()['Total W'].index

    # Get the average poloidal profile from LAMS data where at each "slice" of a
    # bathub curve it is normalized to the max. Then each "slice" is averaged.
    bathtubs = len(otf_df.iloc[otf_df.index<=rad_cutoff*10].index)
    all_tubs_otf = np.zeros((bathtubs, len(otf_df.iloc[0]['Total W'].index)))
    all_tubs_itf = np.zeros((bathtubs, len(itf_df.iloc[0]['Total W'].index)))
    for i in range(0, bathtubs):

        # Normalized bathtub at each radial slice.
        norm_tub_otf = otf_df.iloc[i]['Total W'] / otf_df.iloc[i]['Total W'].max()
        norm_tub_itf = itf_df.iloc[i]['Total W'] / itf_df.iloc[i]['Total W'].max()
        all_tubs_otf[i] = norm_tub_otf
        all_tubs_itf[i] = norm_tub_itf

    # Then the average normalized bathtub.
    avg_tub_otf = all_tubs_otf.mean(axis=0)
    std_tub_otf = all_tubs_otf.std(axis=0)
    avg_tub_itf = all_tubs_itf.mean(axis=0)
    std_tub_itf = all_tubs_itf.std(axis=0)

    # Finally calculate the average peaking ratio.
    mid_itf = int(len(avg_tub_itf) / 2)
    mid_otf = int(len(avg_tub_otf) / 2)
    avg_itf_peak1 = avg_tub_itf[0]  / avg_tub_itf[mid_itf]
    avg_itf_peak2 = avg_tub_itf[-1] / avg_tub_itf[mid_itf]
    avg_otf_peak1 = avg_tub_otf[0]  / avg_tub_otf[mid_otf]
    avg_otf_peak2 = avg_tub_otf[-1] / avg_tub_otf[mid_otf]
    avg_itf_peak = (avg_itf_peak1 + avg_itf_peak2) / 2
    avg_otf_peak = (avg_otf_peak1 + avg_otf_peak2) / 2
    peak_otfitf = avg_otf_peak / avg_itf_peak
    avg_itf_peak_err1 = avg_itf_peak1 * np.sqrt(np.square(std_tub_itf[0] / avg_tub_itf[0]) + np.square(std_tub_itf[mid_itf] / avg_tub_itf[mid_itf]))
    avg_itf_peak_err2 = avg_itf_peak2 * np.sqrt(np.square(std_tub_itf[-1] / avg_tub_itf[-1]) + np.square(std_tub_itf[mid_itf] / avg_tub_itf[mid_itf]))
    avg_otf_peak_err1 = avg_otf_peak1 * np.sqrt(np.square(std_tub_otf[0] / avg_tub_otf[0]) + np.square(std_tub_otf[mid_otf] / avg_tub_otf[mid_otf]))
    avg_otf_peak_err2 = avg_otf_peak2 * np.sqrt(np.square(std_tub_otf[-1] / avg_tub_otf[-1]) + np.square(std_tub_otf[mid_otf] / avg_tub_otf[mid_otf]))
    #avg_itf_peak_err = np.sqrt(np.square(avg_itf_peak_err1) + np.square(avg_itf_peak_err2))
    #avg_otf_peak_err = np.sqrt(np.square(avg_otf_peak_err1) + np.square(avg_otf_peak_err2))
    avg_itf_peak_err = np.sqrt(np.square(avg_itf_peak1 - avg_itf_peak) + np.square(avg_itf_peak2 - avg_itf_peak))
    avg_otf_peak_err = np.sqrt(np.square(avg_otf_peak1 - avg_otf_peak) + np.square(avg_otf_peak2 - avg_otf_peak))
    peak_otfitf_err = peak_otfitf *  np.sqrt(np.square(avg_itf_peak_err / avg_itf_peak) + np.square(avg_otf_peak_err / avg_otf_peak)) / 2
    print("  ITF Peaking: {:.2f} +/- {:.2f}".format(avg_itf_peak, avg_itf_peak_err))
    print("  OTF Peaking: {:.2f} +/- {:.2f}".format(avg_otf_peak, avg_otf_peak_err))
    print("  Peaking OTF/ITF: {:.2f} +/- {:.2f}".format(peak_otfitf, peak_otfitf_err))

    # If there is no RBS data, then we need to also estimate the Total ITF/OTF
    # ratio with LAMS data. This is under the assumption that both have the same
    # calibration curves, which is not necesarilly always the case.
    if p in plist_only_lams:
        itf_totw = itf_df.iloc[itf_df.index<=rad_cutoff*10]["Total W"].sum().sum()
        otf_totw = otf_df.iloc[otf_df.index<=rad_cutoff*10]["Total W"].sum().sum()
        tot_itfotf = itf_totw / otf_totw
        tot_itfotf_err = tot_itfotf * 0.1

    # Sometime there aren't values for the RBS errors so just assign 10%.
    if tot_itfotf_err == 0:
        tot_itfotf_err = tot_itfotf * 0.1

    # Append to lists.
    if p in plist_rev:
        rev_all_tot_itfotf.append(tot_itfotf)
        rev_all_tot_itfotf_err.append(tot_itfotf_err)
        rev_all_peak_otfitf.append(peak_otfitf)
        rev_all_peak_otfitf_err.append(peak_otfitf_err)
    elif p in plist_for:
        for_all_tot_itfotf.append(tot_itfotf)
        for_all_tot_itfotf_err.append(tot_itfotf_err)
        for_all_peak_otfitf.append(peak_otfitf)
        for_all_peak_otfitf_err.append(peak_otfitf_err)
    else:
        print("Error: Assign probe to a Bt direction.")

# Plotting commands.
ms = 14
fig, ax = plt.subplots()
ax.axhline(1.0, linestyle="--", color="k", zorder=1)
ax.axvline(1.0, linestyle="--", color="k", zorder=2)
#ax.errorbar(rev_all_peak_otfitf, rev_all_tot_itfotf, rev_all_tot_itfotf_err, rev_all_peak_otfitf_err, capsize=5, ecolor="k", mec="k", mew=1, fmt='.', ms=ms, color="tab:purple", label=r"Bx$\nabla$B$\uparrow$", zorder=3)
#ax.errorbar(for_all_peak_otfitf, for_all_tot_itfotf, for_all_tot_itfotf_err, for_all_peak_otfitf_err, capsize=5, ecolor="k", mec="k", mew=1, fmt='.', ms=ms, color="tab:red", label=r"Bx$\nabla$B$\downarrow$", zorder=4)
ax.errorbar(rev_all_peak_otfitf, rev_all_tot_itfotf, rev_all_tot_itfotf_err, rev_all_peak_otfitf_err, capsize=5, ecolor="k", mec="k", mew=1, fmt='.', ms=ms, color="tab:purple", label="Unfavorable", zorder=3)
ax.errorbar(for_all_peak_otfitf, for_all_tot_itfotf, for_all_tot_itfotf_err, for_all_peak_otfitf_err, capsize=5, ecolor="k", mec="k", mew=1, fmt='.', ms=ms, color="tab:red", label="Favorable", zorder=4)
ax.plot(lim_otfitf_peak, lim_itfotf_tot, '*', ms=18, color="tab:cyan", mec="k", mew=1, zorder=5, label="3DLIM")
ax.set_xlabel("OTF/ITF Peaking", fontsize=16)
ax.set_ylabel("ITF/OTF Total", fontsize=16)
ax.set_xlim([0, 2])
ax.set_ylim([0, 2])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(fontsize=16, loc="lower right")
fig.tight_layout()
fig.show()
