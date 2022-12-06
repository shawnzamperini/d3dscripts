# A control file to help run through the workflow of BlobbyLP.
from BlobbyLP import BlobbyLP


# The shot and location to the gfile.
#shot = 176487
#time = 3000
#tmin = time
#tmax = 4000
#gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176487/176487_3000.pickle"
#mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176487/lam_xp_length.dat"
#probes = "shelf"

shot = 167195
time = 2500
tmin = 3800
tmax = 5000
gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167195/167195_2500.pickle"
mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_xp_length.dat"
probes = "shelf"

#shot = 174176
#time = 2000
#tmin = 2000
#tmax = 3600
#gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174176/174176_2000.pickle"
#mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174176/lam_xp_length.dat"
#probes = "shelf"

# MAFOT direction we want +1. Since same shapes, can use same parameters here
# as for 190425 and 190427.
#shot = 190423
#time = 3000
#tmin = 2100
#tmax = 3000
#gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190423/190423_3000.pickle"
#mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190423/lam_xp_length.dat"
#probes = "sas"

# Here can be code to optionally load in data from a preset file.
# To-do.

# First we want to fit the TS Te data.
blp = BlobbyLP()
tsdata = blp.calc_te_params(shot, tmin, tmax, tmany=15)
#tsdata = {'tesep': 44.90641076973968,
# 'nesep': 8.514769921852178,
# 'lambda_te': 2.417992167054407,
# 'lambda_ne': 4.105034536977608}

# Then get the flux expansion data.
flux_exp = blp.calc_flux_expansion(shot, time, gfile_pickle_path, target=probes)
#flux_exp = {}

# Then pull the LP data.
lpdata = blp.load_shelf_lp(shot, tmin, tmax, probes=probes)

# Load in the X-point R location.
rxpt = blp.load_xp_r_loc(shot, time)

# Load in the MAFOT connection length data.
mafot = blp.load_mafot(mafot_path)

# Finally run the model. The inputs below are the interpretive model parameters.
# 167195: gengamma, mu_pois=10, vr_mean=20, vr_std=375, gamma_a=1.8, gamma_c=0.80, qtim=1e-7
# 174176: gengamma, mean=20, vr_std=140, gamma_a=1.8, gamma_c=0.75
dist_type = "gengamma"
mu_pois   = 10.0 # poisson only
vr_mean   = 20  # skewnorm, poisson and gamma
vr_std    = 375  # skewnorm, gengamma and gamma only
vr_skew   = 0    # skewnorm only
gamma_a   = 1.8  # gengamma and gamma only
gamma_c   = 0.80  # gengamma only, when 1.0 it's just a gamma dist.
xp_loc    = 5.0  # Not used when MAFOT is loaded
nparts    = 2000
prof_coef = 1.0  # 1.0 = normalize values. Figure this out last.
temin     = 1.0
qtim      = 1e-7
blp.run_blp(tsdata, lpdata, flux_exp, vr_mean, vr_std, vr_skew, xp_loc, nparts,
    prof_coef, dist_type, mu_pois, temin, qtim, gamma_a, gamma_c)

# For using in DIVIMP.
import numpy as np
from scipy.stats import skewnorm
from scipy.interpolate import interp1d

# Get initial distribution with cutoff.
plot_vr = np.linspace(0, vr_mean*100, 10000)
pdf = skewnorm.pdf(plot_vr, a=vr_skew, loc=vr_mean, scale=vr_std)
cutoff = np.where(pdf >= 0.1 * pdf.max())
plot_vr = plot_vr[cutoff]
pdf = pdf[cutoff]

# Resample distribution with smaller number of points.
f = interp1d(plot_vr, pdf)
vrs = np.linspace(plot_vr.min(), plot_vr.max(), 150)
probs = f(vrs)
for i in range(0, len(vrs)):
    print("{:8.2f} {:8.2e}".format(vrs[i], probs[i]))
