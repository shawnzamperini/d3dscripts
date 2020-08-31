import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter


# Flags determining color on plots.
shelf_color = False
bt_color = True
smooth_itfotf = True

# Flags on how to treat the conneciton length data, either apply the shifts in
# probe_shift or cap the connection length at a certain value.
shift = True
cap_conn = 5  # m
min_conn = 0.1

# Don't include data in the wall region.
skip_wall = True

# Use LAMS data "calibrated" with the RBS data.
use_cal_lams = True

# Don't divide by connection length, for comparision.
divide_conn = True

# Smooth the connection length data some to mitigate steps.
smooth_conn = True

# Path prefixes.
probe_root = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/"
shot_root  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/Connection Lengths/"

# Dictionary with the probe and corresponding shot connection length data and
# a dictionary to tell if reverse or forward for ITF/OTF identification. Also
# one saying how many shots it was in for. And where the SP was.
probe_shot = {"A2":"167196", "A3":"167227", "A4":"167229", "A7":"167237", "A8":"167247", "A11":"167266", "A12":"167268", "A15":"167277", "A17":"167279", "A18":"167320", "A19":"167321", "A20":"167322", "A21":"167353", "A33":"167530", "A34":"167534", "A35":"167536"}
probe_side = {"A2":"R", "A3":"R", "A4":"R", "A7":"R", "A8":"R", "A11":"F", "A12":"F", "A15":"F", "A17":"F", "A18":"F", "A19":"F", "A20":"F", "A21":"F", "A33":"F", "A34":"F", "A35":"F"}
probe_nums = {"A2":25, "A3":2, "A4":2, "A7":2, "A8":1, "A11":1, "A12":1, "A15":1, "A17":1, "A18":1, "A19":1, "A20":1, "A21":2, "A33":2, "A34":2, "A35":2}
probe_sp   = {"A2":"S", "A3":"S", "A4":"S", "A7":"S", "A8":"F", "A11":"F", "A12":"F", "A15":"F", "A17":"F", "A18":"S", "A19":"S", "A20":"S", "A21":"S", "A33":"S", "A34":"S", "A35":"S"}
probe_shift = {"A2":-1, "A3":0, "A4":1.5, "A7":-1.0, "A8":-2.5, "A11":0, "A12":0, "A15":-1, "A17":-2, "A18":-3, "A19":-3, "A20":-2, "A21":-2, "A33":-2.5, "A34":-2.5, "A35":-2}  # in cm. neg = probe data shift outwards
probe_cals = {"AD3":2.23E6, "AU3":1.39E6, "AD4":3.27E6, "AU4":1.57E6, "AD8":2.06E6, "AU8":1.79E6, "AD15":2.07E6, "AU15":1.4E6, "AD34":1.45E6, "AU34":1.69E6, "AD35":1.91E6, "AU35":8.7E5}
probe_loc_to_romp = {"A3":(1.1171, 7.5512), "A4":(1.1175, 7.5721), "A8":(1.1241, 8.7533), "A15":(1.0942, 7.1448), "A34":(0.9941, 7.2981), "A35":(0.9982, 7.5731)}

# Locations (mm) along the probe indicating scraped regions.
probe_scrapeD = {"A3":(0,0), "A4":(3.5561,21.08315), "A8":(0,0), "A15":(13.7165,25.4011), "A34":(0,0), "A35":(0,0)}
probe_scrapeU = {"A3":(0,0), "A4":(2.7945,19.305), "A8":(0,0), "A15":(3.30265,16.25705), "A34":(0,0), "A35":(0,0)}

# The location of the wall to baffle jump (i.e. the first baffle point in R-Rsep OMP) from looking at the connection length excel files.
wall_locs  = {"A2":10.20716, "A3":9.96474, "A4":10.69808, "A7":10.1531, "A8":10.13494, "A11":10.13552, "A12":10.30675, "A15":10.09404, "A17":11.99976, "A18":9.4989, "A19":9.49399, "A20":9.80564, "A21":9.01295, "A33":9.3331, "A34":9.26446, "A35":9.3352}
w_dict_itf = {}; w_dict_otf = {}

# Probe to ignore (scraped or such).
p_ignore = ["A4"]
#p_ignore = []

# Probe to plot separately to make sure the W profile and conn line up.
probe_check = ""

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, figsize=(15,5))

for probe, shot in probe_shot.items():

    # Skip ignored probes.
    if probe in p_ignore:
        continue

    print(probe)

    # Load in the data for the probe and corresponding shot.
    shot_itf_df  = pd.read_excel(shot_root + shot + '/' + shot + '.xlsx', sheet_name="MAFOT ITF", skiprows=2)
    shot_otf_df  = pd.read_excel(shot_root + shot + '/' + shot + '.xlsx', sheet_name="MAFOT OTF", skiprows=2)

    # For each probe R-Rsep OMP value find the nearest index in the shot_df.
    # Identify the correct ITF side.
    probe_df = pd.read_excel(probe_root + probe + ".xlsx")
    if probe_side[probe] == "R":
        r_itf = [r for r in probe_df["R-Rsep omp D (cm)"] if not np.isnan(r)]
        r_otf = [r for r in probe_df["R-Rsep omp U (cm)"] if not np.isnan(r)]
        w_itf = [w for w in probe_df["W Areal Density D (1e15 W/cm2)"] if not np.isnan(w)]
        w_otf = [w for w in probe_df["W Areal Density U (1e15 W/cm2)"] if not np.isnan(w)]
    else:
        r_itf = [r for r in probe_df["R-Rsep omp U (cm)"] if not np.isnan(r)]
        r_otf = [r for r in probe_df["R-Rsep omp D (cm)"] if not np.isnan(r)]
        w_itf = [w for w in probe_df["W Areal Density U (1e15 W/cm2)"] if not np.isnan(w)]
        w_otf = [w for w in probe_df["W Areal Density D (1e15 W/cm2)"] if not np.isnan(w)]

    # Overwrite the data with "calibrated" LAMS data (calibrations from just comparing
    # LAMS to RBS).
    if use_cal_lams:

        # Check if LAMS+RBS data was available for this probe.
        probe_u = "AU" + probe.split("A")[1]
        probe_d = "AD" + probe.split("A")[1]
        if probe_u in probe_cals:

            # Get the LAMS data from the Excel sheets. Need to pad the
            # numbers less than 10 with a zero.
            num = int(probe.split("A")[1])
            if num < 10:
                num = "0{}".format(num)
            else:
                num = "{}".format(num)
            lams_pathD  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/AD{}_Map_Analysis.xlsx".format(num)
            lams_pathU  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/AU{}_Map_Analysis.xlsx".format(num)
            lamsD_df = pd.read_excel(lams_pathD, sheet_name="MapData")
            lamsU_df = pd.read_excel(lams_pathU, sheet_name="MapData")

            # Pull out the centerline LAMS data for each side, convert to areal
            # density, and from loc to r-rsep omp. Also replace data in the scraped
            # regions, if specified, with linear fits.
            pol = 2.0
            idxD = lamsD_df["z Location [mm]"] == pol
            lams_locsD = lamsD_df["Axial Location [mm]"][idxD].values
            lams_totWD = lamsD_df["Total W"][idxD].values
            idxU = lamsU_df["z Location [mm]"] == pol
            lams_locsU = lamsU_df["Axial Location [mm]"][idxU].values
            lams_totWU = lamsU_df["Total W"][idxU].values

            # Linear fit across the scraped region.
            if probe_scrapeD[probe] != (0, 0):
                start_scrape_idx = np.where(lams_locsD == probe_scrapeD[probe][0])[0]
                end_scrape_idx   = np.where(lams_locsD == probe_scrapeD[probe][1])[0]
                w_start = lams_totWD[start_scrape_idx]
                w_end = lams_totWD[end_scrape_idx]
                m = (w_end - w_start) / (probe_scrapeD[probe][1] - probe_scrapeD[probe][0])
                scrape_locs_idx = np.logical_and(lams_locsD >= probe_scrapeD[probe][0], lams_locsD <= probe_scrapeD[probe][1])
                w_replace = m * (lams_locsD[scrape_locs_idx] - probe_scrapeD[probe][0]) + w_start
                lams_totWD[scrape_locs_idx] = w_replace
            if probe_scrapeU[probe] != (0, 0):
                dist_start = np.abs(lams_locsU - probe_scrapeU[probe][0])
                start_scrape_idx = np.where(dist_start == dist_start.min())[0]
                dist_end = np.abs(lams_locsU - probe_scrapeU[probe][1])
                end_scrape_idx = np.where(dist_end == dist_end.min())[0]
                #start_scrape_idx = np.where(lams_locsU == probe_scrapeU[probe][0])[0]
                #end_scrape_idx   = np.where(lams_locsU == probe_scrapeU[probe][1])[0]
                w_start = lams_totWU[start_scrape_idx]
                w_end = lams_totWU[end_scrape_idx]
                m = (w_end - w_start) / (probe_scrapeU[probe][1] - probe_scrapeU[probe][0])
                scrape_locs_idx = np.logical_and(lams_locsU >= probe_scrapeU[probe][0], lams_locsU <= probe_scrapeU[probe][1])
                w_replace = m * (lams_locsU[scrape_locs_idx] - probe_scrapeU[probe][0]) + w_start
                lams_totWU[scrape_locs_idx] = w_replace

            # Convert from distance along probe to R-Rsep OMP and counts to areal density.
            lams_locsD = probe_loc_to_romp[probe][0] * lams_locsD / 10 + probe_loc_to_romp[probe][1]
            lams_totWD = lams_totWD / probe_cals[probe_d]
            lams_locsU = probe_loc_to_romp[probe][0] * lams_locsU / 10 + probe_loc_to_romp[probe][1]
            lams_totWU = lams_totWU / probe_cals[probe_u]

            #pol = 2.0
            #idxD = lamsD_df["z Location [mm]"] == pol
            #lams_locsD = probe_loc_to_romp[probe][0] * lamsD_df["Axial Location [mm]"][idxD].values / 10 + probe_loc_to_romp[probe][1] # mm to cm
            #lams_totWD = lamsD_df["Total W"][idxD].values / probe_cals[probe_d]
            #idxU = lamsU_df["z Location [mm]"] == pol
            #lams_locsU = probe_loc_to_romp[probe][0] * lamsU_df["Axial Location [mm]"][idxU].values / 10 + probe_loc_to_romp[probe][1]# mm to cm
            #lams_totWU = lamsU_df["Total W"][idxU].values / probe_cals[probe_u]

            # Assign to correct side, either ITF or OTF.
            if probe_side[probe] == "R":
                r_itf = lams_locsD
                r_otf = lams_locsU
                w_itf = lams_totWD
                w_otf = lams_totWU
            else:
                r_itf = lams_locsU
                r_otf = lams_locsD
                w_itf = lams_totWU
                w_otf = lams_totWD

        # If no RBS+LAMS data, move onto the next one.
        else:
            print("  No RBS+LAMS data available. Skipped.")
            continue


    r_shift = []
    r_conns = []
    conns = []
    w_per_conn_itf = []
    for idx in range(0, len(r_itf)):

        # Get the R values for the connection length data, applying a shift if
        # desired. The shift instead of shifting only one dataset, is split
        # equally between the two. I.e instead of shifting say the RBS data 2 cm
        # outwards, we shift it 1 cm outwards and the conn. data 1 cm inwards.
        if shift:
            #r_conn = shot_itf_df["R-Rsep OMP (cm)"] + probe_shift[probe] / 2
            #r = r_itf[idx] - probe_shift[probe] / 2
            r_conn = shot_itf_df["R-Rsep OMP (cm)"] - probe_shift[probe] / 2
            r = r_itf[idx] + probe_shift[probe] / 2
        else:
            r_conn = shot_itf_df["R-Rsep OMP (cm)"]
            r = r_itf[idx]
        r_shift.append(r)

        # Get nearest index of the shot_df to this R-Rsep OMP value.
        dist = np.abs(r_conn - r)
        close_idx = np.where(dist == np.min(dist))[0][0]

        # Half connection length data, cap at a max value if desired.
        #conn = shot_itf_df["Connection Length (km)"].iloc[close_idx] * 1000 / 2 # km to m
        conn_data = shot_itf_df["Connection Length (km)"]
        if smooth_conn:
            conn_data = savgol_filter(conn_data, 21, 1)
        conn = conn_data[close_idx] * 1000 / 2 # km to m
        if conn > cap_conn:
            conn = cap_conn
        elif conn < min_conn:
            conn = min_conn
        conns.append(conn)

        # If we are in the wall region skip.
        if skip_wall:

            # Get rid of the shift if it was applied since the input data for
            # the wall location is not shifted or anything.
            if shift:
                r_conn_check = r_conn.iloc[close_idx] + probe_shift[probe] / 2
            else:
                r_conn_check = r_conn.iloc[close_idx]

            # If we're in the wall region then scrap the data.
            if r_conn_check > wall_locs[probe]:
                w_per_conn_itf.append(np.nan)
                continue

        # If there is no W collected, just give it a nan.
        if w_itf[idx] == 0:
            w_per_conn_itf.append(np.nan)
            continue

        # Calculate W per connection length per shot.
        if divide_conn:
            w_per_conn_itf.append(w_itf[idx] / conn / probe_nums[probe])
            ytitle = "W Areal Density per Conn. Length per Shot (1e15 W/cm2/m)"
        else:
            w_per_conn_itf.append(w_itf[idx] / probe_nums[probe])
            ytitle = "W Areal Density per per Shot (1e15 W/cm2)"

    # Add to plot, ignoring nans.
    keep_idx = np.where(~np.isnan(w_per_conn_itf))[0]

    # Color according to SP location.
    if shelf_color:
        if probe_sp[probe] == "S":
            color = "tab:purple"
        else:
            color = "tab:red"
        #ax1.plot(np.array(r_itf)[keep_idx], np.array(w_per_conn_itf)[keep_idx], color=color, label=probe)
        ax1.plot(np.array(r_shift)[keep_idx], np.array(w_per_conn_itf)[keep_idx], color=color, label=probe)

    # Color according to Bt direction.
    elif bt_color:
        if probe_side[probe] == "R":
            color = "tab:purple"
        else:
            color = "tab:red"
        #ax1.plot(np.array(r_itf)[keep_idx], np.array(w_per_conn_itf)[keep_idx], color=color, label=probe)
        ax1.plot(np.array(r_shift)[keep_idx], np.array(w_per_conn_itf)[keep_idx], color=color, label=probe)

    # Normal plot.
    else:
        #ax1.plot(np.array(r_itf)[keep_idx], np.array(w_per_conn_itf)[keep_idx], label=probe)
        ax1.plot(np.array(r_shift)[keep_idx], np.array(w_per_conn_itf)[keep_idx], label=probe)

    if probe == probe_check:
        fig_check, ax_check = plt.subplots()
        ax_check2 = ax_check.twinx()
        ax_check.plot(r_shift, w_itf, label=probe)
        plot_conn = []
        for c in shot_itf_df["Connection Length (km)"] / 2 * 1000:
            if c > cap_conn:
                plot_conn.append(cap_conn)
            else:
                plot_conn.append(c)
        ax_check2.semilogy(r_conn, plot_conn, label=shot)
        ax_check.legend()
        ax_check.set_title("#{} - {}".format(shot, probe))
        ax_check.set_ylabel("W Areal Density (1e15 W/cm2)")
        ax_check2.set_ylabel("Connection Length (m)")
        ax_check.set_label("R-Rsep OMP (cm)")
        fig_check.tight_layout()
        fig_check.show()

    # Add to dictionary.
    w_dict_itf[probe] = np.vstack((r_shift, conns, w_per_conn_itf))

    # Repeat for OTF side.
    r_shift = []
    r_conns = []
    conns = []
    w_per_conn_otf = []
    for idx in range(0, len(r_otf)):
        if shift:
            r_conn = shot_otf_df["R-Rsep OMP (cm)"] - probe_shift[probe] / 2
            r = r_otf[idx] + probe_shift[probe] / 2
        else:
            r_conn = shot_otf_df["R-Rsep OMP (cm)"]
            r = r_otf[idx]
        r_shift.append(r)
        dist = np.abs(r_conn - r)
        close_idx = np.where(dist == np.min(dist))[0][0]
        #conn = shot_otf_df["Connection Length (km)"].iloc[close_idx] * 1000 / 2 # km to m
        conn_data = shot_otf_df["Connection Length (km)"]
        if smooth_conn:
            conn_data = savgol_filter(conn_data, 21, 1)
        conn = conn_data[close_idx] * 1000 / 2 # km to m
        if conn > cap_conn:
            conn = cap_conn
        conns.append(conn)
        if w_otf[idx] == 0:
            w_per_conn_otf.append(np.nan)
            continue

        if divide_conn:
            w_per_conn_otf.append(w_otf[idx] / conn / probe_nums[probe])
            ytitle = "W Areal Density per Conn. Length per Shot (1e15 W/cm2/m)"
        else:
            w_per_conn_otf.append(w_otf[idx] / probe_nums[probe])
            ytitle = "W Areal Density per per Shot (1e15 W/cm2)"

        #w_per_conn_otf.append(w_otf[idx] / conn / probe_nums[probe])
    keep_idx = np.where(~np.isnan(w_per_conn_otf))[0]
    if shelf_color:
        if probe_sp[probe] == "S":
            color = "tab:purple"
        else:
            color = "tab:red"
        ax2.plot(np.array(r_shift)[keep_idx], np.array(w_per_conn_otf)[keep_idx], color=color, label=probe)
    elif bt_color:
        if probe_side[probe] == "R":
            color = "tab:purple"
        else:
            color = "tab:red"
        ax2.plot(np.array(r_shift)[keep_idx], np.array(w_per_conn_otf)[keep_idx], color=color, label=probe)
    else:
        ax2.plot(np.array(r_shift)[keep_idx], np.array(w_per_conn_otf)[keep_idx], label=probe)
    if probe == probe_check:
        fig_check, ax_check = plt.subplots()
        ax_check2 = ax_check.twinx()
        ax_check.plot(r_shift, w_otf, label=probe)
        ax_check2.semilogy(r_conn, shot_otf_df["Connection Length (km)"], label=shot)
        ax_check.legend()
        ax_check.set_title("#{} - {}".format(shot, probe))
        ax_check.set_ylabel("W Areal Density (1e15 W/cm2)")
        ax_check2.set_ylabel("Connection Length (km)")
        ax_check.set_label("R-Rsep OMP (cm)")
        fig_check.tight_layout()
        #fig_check.show()
    w_dict_otf[probe] = np.vstack((r_shift, conns, w_per_conn_otf))

"""
    # Repeat for OTF side.
    w_per_conn_otf = []
    for idx in range(0, len(r_otf)):
        if shift:
            r_conn = shot_otf_df["R-Rsep OMP (cm)"] + probe_shift[probe] / 2
            r = r_otf[idx] - probe_shift[probe] / 2
        else:
            r_conn = shot_otf_df["R-Rsep OMP (cm)"]
            r = r_otf[idx]
        dist = np.abs(r_conn - r)
        close_idx = np.where(dist == np.min(dist))[0][0]
        conn = shot_otf_df["Connection Length (km)"].iloc[close_idx] * 1000 / 2  # km to m
        w_per_conn_otf.append(w_otf[idx] / conn / probe_nums[probe])
    keep_idx = np.where(~np.isnan(w_per_conn_otf))[0]
    ax2.plot(np.array(r_otf)[keep_idx], np.array(w_per_conn_otf)[keep_idx], label=probe)
    w_dict_otf[probe] = np.vstack((r_otf, w_per_conn_otf))
"""

# Loop through our output dictionaries and use spline fits to calculate the
# ITF/OTF ratios for each probe to plot on ax3.
for probe, data in w_dict_itf.items():

    # Bounds for interpolation.
    rmin = np.max((w_dict_itf[probe][0].min(), w_dict_otf[probe][0].min()))
    rmax = np.min((w_dict_itf[probe][0].max(), w_dict_otf[probe][0].max()))
    r_int = np.linspace(rmin, rmax, 1000)

    # Interpolation functions for each side.
    itf_int = interp1d(w_dict_itf[probe][0], w_dict_itf[probe][2])
    otf_int = interp1d(w_dict_otf[probe][0], w_dict_otf[probe][2])
    itfotf = itf_int(r_int) / otf_int(r_int)

    if smooth_itfotf:
        itfotf = savgol_filter(itfotf, 51, 3)

    # Plotting commands.

    if bt_color:
        if probe_side[probe] == "R":
            color = "tab:purple"
        else:
            color = "tab:red"
        ax3.plot(r_int, itfotf, label=probe, color=color)
    else:
        ax3.plot(r_int, itfotf, label=probe)
    #ax3.semilogy(r_int, itfotf, label=probe)

ax1.set_ylabel(ytitle)
ax1.set_xlabel("R-Rsep OMP (cm)")
ax2.set_xlabel("R-Rsep OMP (cm)")
ax3.set_xlabel("R-Rsep OMP (cm)")
ax3.set_ylabel("ITF/OTF")
#ax1.set_ylim([0, 0.04])
ax2.set_ylim([0, 0.04])
ax3.set_xlim([None, 12])
ax3.set_ylim([0, None])
ax3.axhline(1, color="k", linestyle="--")
#ax1.legend()
ax2.legend()
fig.tight_layout()
fig.show()
