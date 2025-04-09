import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
import matplotlib.ticker as ticker


# Figure of merit. This is what is used to compare the different scenarios
# to each other. Options are: avg_conc, avg_nw
fom = "avg_conc"

# Path to directory where files are stored
root = "/fusion/projects/codes/oedge/zamperinis/results"

# Choose the ring just inside the separatrix for analysis
ring = 15

# Let's assume the divertors have a high W prompt redeposition fraction of about 99%,
# and everywhere else has a 50% chance.
fprompt_w = {"all":0.0, "sh":0.99, "id":0.99, "ub":0.50, "iw":0.50, "ow":0.50, "lb":0.50}
#fprompt_w = {"all":0.0, "sh":0.99, "id":0.99, "ub":0.00, "iw":0.00, "ow":0.00, "lb":0.00}

# Just half the W number for Fe
fprompt_fe = {k:v/2 for k,v in fprompt_w.items()}

# Create figure for the plotting below
fontsize = 16
fig, ax = plt.subplots()

# The parameters that were scanned
dperps = [3.098, 2.745, 0.165, 3.169, 4.932, 1.959, 4.663, 0.544, 1.065, 1.514, 4.965, 4.872, 0.774, 3.235, 4.632, 3.682, 0.576, 1.990, 3.747, 0.688, 3.949, 4.285, 3.763, 4.893, 1.131, 4.231, 2.197, 2.344, 0.846, 4.469, 4.910, 3.191, 2.007, 0.328, 3.910, 2.863, 1.388, 0.440, 2.593, 3.309, 4.882, 1.696, 3.468, 3.160, 4.345, 0.498, 3.611, 1.132, 3.643, 0.260, 0.851, 0.286, 4.777, 3.226, 4.065, 0.440, 2.848, 0.489, 3.471, 2.578, 0.103, 4.154, 0.768, 3.622, 4.592, 2.799, 2.716, 1.198, 3.903, 2.957, 4.026, 3.158, 3.727, 0.501, 2.621, 4.008, 4.689, 4.445, 1.981, 4.851, 4.073, 3.174, 4.139, 1.290, 0.788, 0.439, 3.655, 4.455, 1.152, 4.842, 1.576, 0.941, 3.938, 0.897, 4.591, 0.995, 0.746, 2.619, 3.154, 3.690]
ne_charges = [4, 8, 7, 7, 5, 7, 3, 6, 6, 5, 2, 5, 4, 3, 8, 9, 8, 8, 6, 2, 2, 2, 3, 9, 4, 8, 2, 2, 4, 5, 7, 2, 8, 2, 7, 7, 2, 4, 4, 2, 4, 3, 2, 9, 2, 2, 5, 3, 7, 3, 7, 5, 3, 8, 4, 2, 9, 9, 5, 6, 6, 9, 4, 7, 9, 6, 7, 8, 3, 9, 4, 4, 4, 3, 8, 7, 8, 8, 3, 8, 9, 2, 8, 2, 5, 9, 2, 7, 5, 4, 8, 6, 6, 8, 5, 9, 3, 4, 2, 7]
ne_fluxes = [0.052, 0.025, 0.035, 0.059, 0.054, 0.044, 0.082, 0.085, 0.006, 0.035, 0.100, 0.010, 0.053, 0.044, 0.080, 0.027, 0.067, 0.053, 0.050, 0.036, 0.036, 0.046, 0.039, 0.039, 0.059, 0.032, 0.095, 0.015, 0.054, 0.092, 0.055, 0.082, 0.083, 0.047, 0.055, 0.031, 0.012, 0.009, 0.091, 0.055, 0.028, 0.100, 0.097, 0.040, 0.060, 0.075, 0.023, 0.028, 0.059, 0.082, 0.006, 0.040, 0.036, 0.057, 0.087, 0.012, 0.040, 0.088, 0.098, 0.021, 0.056, 0.053, 0.012, 0.082, 0.054, 0.089, 0.049, 0.060, 0.062, 0.032, 0.056, 0.084, 0.065, 0.025, 0.019, 0.021, 0.016, 0.013, 0.054, 0.082, 0.070, 0.081, 0.044, 0.096, 0.074, 0.016, 0.090, 0.059, 0.065, 0.065, 0.097, 0.036, 0.078, 0.046, 0.042, 0.025, 0.051, 0.063, 0.077, 0.029]
wall_tes = [6.463, 6.144, 4.056, 9.486, 7.510, 2.143, 8.307, 9.969, 9.983, 4.563, 5.143, 5.467, 9.094, 6.464, 4.083, 6.577, 8.594, 7.436, 2.875, 8.249, 1.022, 7.952, 5.742, 5.539, 3.894, 2.760, 1.057, 1.455, 6.794, 4.243, 1.001, 8.149, 8.091, 1.217, 8.500, 5.004, 1.846, 9.438, 3.976, 1.852, 3.899, 7.828, 9.845, 5.489, 5.287, 9.320, 2.102, 4.305, 2.908, 3.358, 2.362, 8.786, 8.685, 7.122, 2.042, 9.576, 4.264, 2.809, 2.321, 5.893, 8.132, 6.327, 1.381, 9.911, 8.369, 6.640, 4.480, 2.999, 4.433, 3.707, 7.064, 2.570, 8.686, 4.260, 4.540, 8.724, 3.106, 4.405, 9.789, 6.546, 7.609, 5.961, 3.606, 1.352, 8.419, 4.236, 6.324, 5.655, 3.145, 9.622, 9.333, 1.480, 5.782, 2.052, 7.440, 3.836, 8.712, 8.237, 5.187, 2.251]
machs = [-0.210, -0.381, -0.416, -0.338, -0.001, -0.050, -0.221, -0.243, -0.118, -0.203, -0.003, -0.071, -0.418, -0.067, -0.430, -0.397, -0.205, -0.146, -0.180, -0.230, -0.372, -0.433, -0.405, -0.139, -0.086, -0.423, -0.125, -0.196, -0.497, -0.319, -0.371, -0.059, -0.332, -0.016, -0.421, -0.373, -0.114, -0.067, -0.280, -0.417, -0.255, -0.199, -0.064, -0.339, -0.248, -0.487, -0.186, -0.116, -0.392, -0.152, -0.121, -0.374, -0.087, -0.285, -0.286, -0.147, -0.015, -0.384, -0.243, -0.377, -0.410, -0.240, -0.303, -0.128, -0.371, -0.004, -0.085, -0.278, -0.462, -0.056, -0.329, -0.041, -0.313, -0.052, -0.147, -0.218, -0.075, -0.130, -0.470, -0.475, -0.218, -0.463, -0.408, -0.455, -0.085, -0.023, -0.189, -0.064, -0.326, -0.226, -0.390, -0.478, -0.115, -0.322, -0.258, -0.450, -0.185, -0.109, -0.225, -0.416]

# Zip together so we have a record of each case's setting
xvalues = np.array(list(zip(dperps, ne_charges, ne_fluxes, wall_tes, machs)))
#xvalues = xvalues[:3]

# Make Mach values positive, arbitrary
xvalues[:,4] *= -1

# Normalize the xvalues between 0-1
#for i in range(0, xvalues.shape[1]):
#   xvalues[:,i] = (xvalues[:,i] - xvalues[:,i].min()) / \
#       (xvalues[:,i].max() - xvalues[:,i].min())

# Normalizing data
scaler = StandardScaler()
xvalues_norm = scaler.fit_transform(xvalues)

# Some quick 2D plots that are used in the presentation
ncpath = "{}/d3d-w-wall-param-scan-v5-w-all-1.nc".format(root)
op = oedge_plots.OedgePlots(ncpath)
ktebs = op.read_data_2d("KTEBS")
ones = np.ones(ktebs.shape)
op.plot_contour_polygon("KTEBS", own_data=ones, vmin=0, vmax=2, cmap="coolwarm")
#op.plot_contour_polygon("KTEBS", normtype="log", cbar_label="Te (eV)", 
#    cmap="inferno", vmin=1)
#op.plot_contour_polygon("DDLIMS", normtype="log", cbar_label="W Density (m-3)",
#    charge="all", cmap="inferno", scaling=op.absfac, vmin=1e16, vmax=1e19)

# OpenADAS interface for calculating line radiation. First import it.
import sys
sys.path.append("/home/zamperinis/d3dscripts/Utilities")
import openadas

# Create OpenADAS object 
oa = openadas.OpenADAS()

# DataFrames for line radiation (plt) and bremstrahlung (prb)
w_plt = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/plt50/plt50_w.dat")
w_prb = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/prb50/prb50_w.dat")
fe_plt = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/plt89/plt89_fe.dat")
fe_prb = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/prb89/prb89_fe.dat")

# Radiation values caluclated for a typical core plasma
core_te = 5e3
core_ne = 1e20 # Is this typical?
core_rad = 0.20  # Radius wherein these plasma parameters exist
core_vol = core_rad**2 * np.pi

# Get the respective powers. I can only assume these are W m3 (they've been internally
# multiplied within the OpenADAS class by 1e-6 to do cm3 --> m3). 
w_plt_rates = np.array([oa.get_rate_coef(w_plt, core_te, core_ne, charge=z+1) for z in range(74)])
w_prb_rates = np.array([oa.get_rate_coef(w_prb, core_te, core_ne, charge=z+1) for z in range(74)])
fe_plt_rates = np.array([oa.get_rate_coef(fe_plt, core_te, core_ne, charge=z+1) for z in range(26)])
fe_prb_rates = np.array([oa.get_rate_coef(fe_prb, core_te, core_ne, charge=z+1) for z in range(26)])

# DataFrames for ionization (scd) and recombination (acd) rates
w_scd = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/scd50/scd50_w.dat")
w_acd = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/acd50/acd50_w.dat")
fe_scd = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/scd89/scd89_fe.dat")
fe_acd = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/acd89/acd89_fe.dat")

# Get the fractional abundances of each charge state for the core parameters
w_frac_abund = oa.coronal_equil_frac_abundance(74, w_scd, w_acd, core_te, core_ne)
fe_frac_abund = oa.coronal_equil_frac_abundance(26, w_scd, w_acd, core_te, core_ne)

fom_vals_dict = {}
fom_stds_dict = {}
#sput_regions = ["sh", "id", "ub", "ow", "lb", "iw"]
sput_regions = ["all", "ub", "ow", "lb", "sh", "iw", "id"]
for sput_region in sput_regions:

    # Lists to hold our figures of merit as we extract them from each case
    fom_vals = []
    fom_stds = []
    w_conc = []
    fe_conc = []
    w_prad = []
    fe_prad = []

    for i in range(0, xvalues.shape[0]):

        # fe-comp is to directly compare against the w cases, it is the same
        # surfaces coated in each, respectively
        for imp in ["w", "fe-comp"]:

            # Load case
            ncpath = ("{}/d3d-w-wall-param-scan-v5-{}-{}-{}.nc"
                .format(root, imp, sput_region, i+1))
            print(ncpath)
            try:
                op = oedge_plots.OedgePlots(ncpath)

            # This can happen if no particles were launched
            except:
                print("  Does not exist! Assigning zeros.")
                fom_vals.append(0)
                fom_stds.append(0)
                if imp == "w":
                    w_conc.append(0)
                    w_prad.append(0)
                else:
                    fe_conc.append(0)
                    fe_prad.append(0)
                continue


            # Load needed density to derive figure of merit at desired ring
            s, nz_s = op.along_ring(ring, "DDLIMS", charge="all", 
                scaling=op.absfac, plot_it=False)
            s, ne_s = op.along_ring(ring, "KNBS", plot_it=False)

            # Multiply by the prompt redeposition fraction as a rough guestimate.
            if imp == "w":
                nz_s *= (1.0 - fprompt_w[sput_region])
            else:
                nz_s *= (1.0 - fprompt_fe[sput_region])

            cz_s = nz_s / (nz_s + ne_s)

            # Now that we have the average concentration, we assume it is constant 
            # throughout the core. Then using the fractional abundances already
            # calculated for typical core conditions, we calculated the radiated
            # power density.
            # To calculate W/m3 it's rate[Z] * ne * nz[Z] = rate[Z] * ne * (ne * w_conc * frac[Z])
            if imp == "w":
                plt_rates = w_plt_rates
                prb_rates = w_prb_rates
                frac_abund = w_frac_abund
            else:
                plt_rates = fe_plt_rates
                prb_rates = fe_prb_rates
                frac_abund = fe_frac_abund

            z_plt_rad = plt_rates * core_ne ** 2 * cz_s.mean() * frac_abund[1:]  # [:1] because skip neutrals
            z_prb_rad = prb_rates * core_ne ** 2 * cz_s.mean() * frac_abund[1:]  # [:1] because skip neutrals

            # Total is just sum of these two contributions across each charge state
            if imp == "w":
                w_prad.append(np.sum(z_plt_rad + z_prb_rad))

                # Figure of merit: Average W concentration. Right now this is a remnant of the previous
                # script. We will probably remove...
                fom_vals.append(cz_s.mean())
                fom_stds.append(cz_s.std())
                w_conc.append(cz_s.mean())
            else:
                fe_prad.append(np.sum(z_plt_rad + z_prb_rad))
                fe_conc.append(cz_s.mean())


    fom_vals_dict[sput_region] = np.array(fom_vals)
    fom_stds_dict[sput_region] = np.array(fom_stds)
    fom_vals_dict[sput_region + "_w_prad"] = np.array(w_prad)
    fom_vals_dict[sput_region + "_fe_prad"] = np.array(fe_prad)
    fom_vals_dict[sput_region + "_w_conc"] = np.array(w_conc)
    fom_vals_dict[sput_region + "_fe_conc"] = np.array(fe_conc)

#ax.set_yscale("log")
ax.set_xlabel("Distance from inner target (m)", fontsize=fontsize)

# Format y axis based off figure of merit.
if fom == "avg_conc":
    ax.set_ylabel("W Concentration", fontsize=fontsize)
    #ax.set_ylim([0, 1e-5])
elif fom == "avg_nw":
    ax.set_ylabel("W Density (m-3)", fontsize=fontsize)

#ax.set_yscale("log")
#ax.set_ylim([5e-6, 3e-2])
fig.tight_layout()
fig.show()

# Map the abbreviations to the plot friendly words
reg_words = {"ib":"Inner\nBaffle", "ub":"Upper\nBaffle", "sh":"Shelf", 
    "ow":"Outer\nWall", "iw":"Inner\nWall", "id":"Inner\nDivertor", 
    "lb":"Lower\nBaffle", "all":"Everywhere"}

# For each region, find the best and worst performers among the set of parameters,
# based off total Prad at core conditions
prad_dict = {}
bar_regions = []
min_prads = {"W":[], "Fe":[]}
avg_prads = {"W":[], "Fe":[]}
max_prads = {"W":[], "Fe":[]}
for sput_region in sput_regions:

    # This doesn't make sense for the "all" cases
    #if sput_region == "all": continue

    # Shelf will already be coated in W, and it throws the plot
    # off anyways so we can just ignore
    #if sput_region == "shelf": continue

    # Find lowest and highest Prad cases
    w_prad = fom_vals_dict[sput_region + "_w_prad"]
    fe_prad = fom_vals_dict[sput_region + "_fe_prad"]
    tot_prad = w_prad + fe_prad
    min_prad_idx = np.argmin(tot_prad)
    max_prad_idx = np.argmax(tot_prad)

    # Save in dictionary that will be used for plotting
    prad_dict[sput_region + "_min_prad"] = tot_prad[min_prad_idx]
    prad_dict[sput_region + "_min_prad_w"] = w_prad[min_prad_idx]
    prad_dict[sput_region + "_min_prad_fe"] = fe_prad[min_prad_idx]
    prad_dict[sput_region + "_max_prad"] = tot_prad[max_prad_idx]
    prad_dict[sput_region + "_max_prad_w"] = w_prad[max_prad_idx]
    prad_dict[sput_region + "_max_prad_fe"] = fe_prad[max_prad_idx]

    # We want really the ratio of W/Fe Prad so we can see when the
    # impact on the core will be worse
    prad_dict[sput_region + "_min_w/fe_Prad"] = w_prad[min_prad_idx] / fe_prad[min_prad_idx]
    prad_dict[sput_region + "_max_w/fe_Prad"] = w_prad[max_prad_idx] / fe_prad[max_prad_idx]

    # Assemble into arrays easy for plotting
    bar_regions.append(sput_region)
    min_prads["W"].append(w_prad[min_prad_idx])
    min_prads["Fe"].append(fe_prad[min_prad_idx])
    avg_prads["W"].append(np.mean(w_prad))
    avg_prads["Fe"].append(np.mean(fe_prad))
    max_prads["W"].append(w_prad[max_prad_idx])
    max_prads["Fe"].append(fe_prad[max_prad_idx])

# Convert to numpy arrays and from MW/m3 to MW
for key, vals in min_prads.items():
    min_prads[key] = np.array(vals) / 1e6 * core_vol
for key, vals in avg_prads.items():
    avg_prads[key] = np.array(vals) / 1e6 * core_vol
for key, vals in max_prads.items():
    max_prads[key] = np.array(vals) / 1e6 * core_vol

# Ratio of W/Fe 
avg_prad_wfe = avg_prads["W"] / avg_prads["Fe"]

# List of the regions as real words
bar_regions_words = [reg_words[k] for k in bar_regions]

# Option to only plot select bars (useful for presentations when
# stepping through bar chart). Order corresponds to sput_regions
#              ["all", "sh", "id", "ub", "ow", "lb", "iw"]
no_plot_bars = [True, False, False, False, False, False, False]

# Make bar chart
width = 0.2
fontsize = 20
colors = ["lightskyblue", "hotpink"] 
y = np.arange(1, len(bar_regions)+1)
fig, ax1 = plt.subplots(figsize=(8, 7))

# Same bar plot, but attempt to split the x axis over multiple values
fig2, axs = plt.subplots(1, 3, sharey=True, figsize=(8, 6))
axs = axs.flatten()

# Put them all in a list, since we will be using the same plotting
# commands on all of them.
all_axs = list(axs)
all_axs.append(ax1)

count = 0
bottom = np.zeros(len(bar_regions))
for key, vals in min_prads.items():
    tmp = vals
    tmp[no_plot_bars] = 0.0
    for ax in all_axs:
        p = ax.barh(y-0.25, tmp, width, left=bottom, color=colors[count], 
            alpha=0.3, edgecolor="k")
    bottom += tmp
    count += 1

count = 0
bottom = np.zeros(len(bar_regions))
for key, vals in avg_prads.items():
    tmp = vals
    tmp[no_plot_bars] = 0.0
    for ax in all_axs:
        p = ax.barh(y, tmp, width, label=key, left=bottom, color=colors[count], 
            tick_label=bar_regions_words, edgecolor="k")
    bottom += tmp
    count += 1

count = 0
bottom = np.zeros(len(bar_regions))
for key, vals in max_prads.items():
    tmp = vals
    tmp[no_plot_bars] = 0.0
    for ax in all_axs:
        p = ax.barh(y+0.25, tmp, width, left=bottom, color=colors[count], 
            alpha=0.3, edgecolor="k")
    bottom += tmp
    count += 1

# Commands for the usual plot
lims = [0, 5]
#lims = [1, 300]
#ax1.set_xscale("log")
ax1.set_xlim(lims)
if core_ne != 1e20:
    print("WRONG XLABEL, FIX NE ASSUMING 1E20")
ne_str = r"$\mathdefault{10^{20}\ m^{-3}}$"
xlabel = "Prad @ {} keV ".format(int(core_te/1e3)) + ne_str + " (MW)"
ax1.set_xlabel(xlabel, fontsize=fontsize)
ax1.tick_params(axis="both", which="major", labelsize=fontsize-2)
ax1.legend(fontsize=fontsize)
fig.tight_layout()

# Commands for split x axis
axs[0].set_xlim([0, 7])
axs[1].set_xlim([12 ,19])
axs[2].set_xlim([43, 50])
axs[0].spines.right.set_visible(False)
axs[1].spines.left.set_visible(False)
axs[1].spines.right.set_visible(False)
axs[2].spines.left.set_visible(False)
axs[1].tick_params(left=False)
axs[2].tick_params(left=False)
d = .0  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-d, -1), (d, 1)], markersize=12, linestyle="none", 
    color='k', mec='k', mew=1, clip_on=False)
axs[0].plot([1, 1], [0, 1], transform=axs[0].transAxes, **kwargs)
axs[1].plot([0, 0], [0, 1], transform=axs[1].transAxes, **kwargs)
axs[1].plot([1, 1], [0, 1], transform=axs[1].transAxes, **kwargs)
axs[2].plot([0, 0], [0, 1], transform=axs[2].transAxes, **kwargs)
axs[0].xaxis.set_major_locator(ticker.MultipleLocator(3))
axs[1].xaxis.set_major_locator(ticker.MultipleLocator(3))
axs[2].xaxis.set_major_locator(ticker.MultipleLocator(3))
axs[1].set_xlabel(xlabel, fontsize=fontsize)
axs[0].tick_params(axis="both", which="major", labelsize=fontsize-2)
axs[1].tick_params(axis="both", which="major", labelsize=fontsize-2)
axs[2].tick_params(axis="both", which="major", labelsize=fontsize-2)
axs[2].legend(fontsize=fontsize)
fig2.tight_layout()
fig2.subplots_adjust(wspace=0.1)

# Horizontal lines (caps) at where bars are split. Buffer just prevents
# half the caps form being cut off since they straddle the edge of the plots
cap_width = 0.14
buf = 0.02
axs[0].plot([7-buf, 7-buf], [5.25-cap_width, 5.25+cap_width], lw=2, color="k", alpha=0.3)
axs[1].plot([12, 12], [5.25-cap_width, 5.25+cap_width], lw=2, color="k", alpha=0.3)
"""
axs[0].plot([7-buf, 7-buf], [2.25-cap_width, 2.25+cap_width], lw=2, color="k", alpha=0.3)
axs[1].plot([12, 12], [2.25-cap_width, 2.25+cap_width], lw=2, color="k", alpha=0.3)
axs[0].plot([7-buf, 7-buf], [1.25-cap_width, 1.25+cap_width], lw=2, color="k", alpha=0.3)
axs[1].plot([12, 12], [1.25-cap_width, 1.25+cap_width], lw=2, color="k", alpha=0.3)
axs[1].plot([19-buf*2, 19-buf*2], [1.25-cap_width, 1.25+cap_width], lw=2, color="k", alpha=0.3)
axs[2].plot([43, 43], [1.25-cap_width, 1.25+cap_width], lw=2, color="k", alpha=0.3)
axs[0].plot([7-buf, 7-buf], [1.00-cap_width, 1.00+cap_width], lw=2, color="k")
axs[1].plot([12, 12], [1.00-cap_width, 1.00+cap_width], lw=2, color="k")
"""

plt.show()

# Plot of average ratio of W/Fe
fig, ax1 = plt.subplots()
ax1.barh(y, avg_prad_wfe, 0.75, tick_label=bar_regions_words, edgecolor="k")
ax.set_xlabel("W/Fe")
fig.tight_layout()
plt.show()

