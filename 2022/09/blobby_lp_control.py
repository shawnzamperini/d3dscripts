# A control file to help run through the workflow of BlobbyLP.
from BlobbyLP import BlobbyLP


# The shot and location to the gfile.
shot = 176487
time = 3000
tmin = time
tmax = 4000
gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176487/176487_3000.pickle"
mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176487/lam_xp_length.dat"

# Here can be code to optionally load in data from a preset file.
# To-do.

# First we want to fit the TS Te data.
blp = BlobbyLP()
#tsdata = blp.calc_te_params(shot, tmin, tmax, tmany=15)
tsdata = {'tesep': 44.90641076973968,
 'nesep': 8.514769921852178,
 'lambda_te': 2.417992167054407,
 'lambda_ne': 4.105034536977608}

# Then get the flux expansion data.
flux_exp = blp.calc_flux_expansion(shot, time, gfile_pickle_path)
#flux_exp = {}

# Then pull the LP data.
lpdata = blp.load_shelf_lp(shot, tmin, tmax)

# Load in the X-point R location.
rxpt = blp.load_xp_r_loc(shot, time)

# Load in the MAFOT connection length data.
mafot = blp.load_mafot(mafot_path)

# Finally run the model. The inputs below are the interpretive model parameters.
dist_type = "skewnorm"
mu_pois = 1.0  # poisson only
vr_mean = 100
vr_std = 750  # skewnorm only
vr_skew = 20   # skewnorm only
xp_loc = 5.0   # Not used when MAFOT is loaded
nparts = 10000
prof_coef = 1.0  # 1.0 = normalize values.
temin = 1.0
qtim = 5e-8
blp.run_blp(tsdata, lpdata, flux_exp, vr_mean, vr_std, vr_skew, xp_loc, nparts,
    prof_coef, dist_type, mu_pois, temin)
