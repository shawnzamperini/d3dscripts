import oedge_plots
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Some constants.
a = 0.585  # m
k = 1.720
sep_ring = 34
last_ring = 70
e = 1.609e-19
boltz = 8.617*1e-5  # eV/K
bt_mult = 1

# Load the 167247 plasma background.
ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-bkg-011f.nc'
op = oedge_plots.OedgePlots(ncpath)

# Load related arrays and constanst from case.
Te = op.read_data_2d("KTEBS")
ne = op.read_data_2d("KNBS")
R = op.read_data_2d("RS")
Z = op.read_data_2d("ZS")
Bt = op.read_data_2d("BTS") * bt_mult
Bp = op.read_data_2d("BRATIO") * Bt
B = np.sqrt(np.square(Bt) + np.square(Bp)) * bt_mult
Er = op.read_data_2d("E_RAD")
R0 = float(op.nc.variables["R0"][:].data)
Z0 = float(op.nc.variables["Z0"][:].data)
pi = ne * Te

# Calculate the radial ion pressure gradient. This is tricky, because Ti is not
# very exponentially decreasing or anything like that and can thus have super
# steep and thus really unrealistic gradients (not to mention it'd be really
# to calculate!). As an approximation, we will just take the Ti, ne values
# for each "row" of knots at the separatrix and at the last ring and draw an
# exponential curve between them and use that for the gradient so it's smooth.
# This probably isn't good below the X-point, so the analysis will be restricted
# to above it, which is fine.
Ti_2d = op.nc.variables["KTIBS"][:].data
ne_2d = op.nc.variables["KNBS"][:].data
R_2d = op.nc.variables["RS"][:].data
Z_2d = op.nc.variables["ZS"][:].data
pi_2d = np.zeros(Ti_2d.shape)
pi_fit_2d = np.zeros(Ti_2d.shape)
pi_grad_2d = np.zeros(Ti_2d.shape)

def exp_fit(x, a, b):
    return a * np.exp(-b * x)

# Loop through each knot index, getting those knots for every ring (here dubbed
# a row of knots). The radial pressure gradient we care about is just in the SOL,
# which we've indexed with sep_ring anf far_ring.
count = 0
for row_knots_idx in range(0, Ti_2d.shape[1]):
    row_knots_Ti = Ti_2d[:, row_knots_idx]
    row_knots_ne = ne_2d[:, row_knots_idx]
    row_knots_R  = R_2d[:, row_knots_idx]
    row_knots_Z  = Z_2d[:, row_knots_idx]

    sep_R  = row_knots_R[sep_ring];  far_R  = row_knots_R[last_ring]
    sep_Z  = row_knots_Z[sep_ring];  far_Z  = row_knots_Z[last_ring]
    row_knots_dist = np.sqrt(np.square(sep_R - row_knots_R) + np.square(sep_Z - row_knots_Z))

    # Calculate ~radial distance from separatrix.
    dist = np.sqrt(np.square(sep_R - far_R) + np.square(sep_Z - far_Z))

    # Calculate ion pressure.
    row_knots_pi = row_knots_ne * row_knots_Ti
    sep_pi = row_knots_pi[sep_ring]; far_pi = row_knots_pi[last_ring]
    pi_2d[:, row_knots_idx] = row_knots_pi

    # Exponential fit between two points for pi.
    popt_pi, pcov = curve_fit(exp_fit, (0, dist), (sep_pi*1e-20, far_pi*1e-20))

    # Fill in the pi gradient array at this "row" of knots.
    row_knots_pi_fit = exp_fit(row_knots_dist, *popt_pi) * 1e20
    row_knots_pi_grad = -popt_pi[1] * exp_fit(row_knots_dist, *popt_pi) * 1e20
    pi_grad_2d[:, row_knots_idx] = row_knots_pi_grad
    pi_fit_2d[:, row_knots_idx] = row_knots_pi_fit

    # Print out to see how things are looking.
    #print("{}: {}".format(count, popt_pi))
    count += 1

# Convert pi_grad array to a 1D array for plotting.
pi_grad = np.zeros(op.num_cells)
pi_fit  = np.zeros(op.num_cells)
count = 0
for ir in range(op.nrs):
    for ik in range(op.nks[ir]):
        if op.area[ir, ik] != 0.0:
            pi_grad[count] = pi_grad_2d[ir][ik]
            pi_fit[count]  = pi_fit_2d[ir][ik]
            count = count + 1

# Calculate the estimated parallel flow speed in sections. First calculate alpha.
alpha = np.arccos(np.abs(R-R0) / np.sqrt(np.square(R-R0) + np.square(Z-Z0)))

# Trig stuff to get real angle.
quad_1 = np.logical_and((R - R0) < 0, (Z - Z0) > 0)
quad_4 = np.logical_and((R - R0) < 0, (Z - Z0) < 0)
quad_3 = np.logical_and((R - R0) > 0, (Z - Z0) < 0)
alpha[quad_1] = np.pi - alpha[quad_1]
alpha[quad_4] = alpha[quad_4] + np.pi
alpha[quad_3] = 2 * np.pi - alpha[quad_3]

# Want to subtract the pf_part on the outside, add on inside.
mult = np.full(alpha.shape, 1)
mult[R > R0] = -1

# Finally get tan^2(alpha) and calculate vcomb.
tan2 = np.square(np.tan(alpha))
exb_part = Er / Bp
pf_part = 2 * (a/R) * (B/Bp) * (pi_grad/(e*ne*B)) * 1 / np.square(1 + tan2/np.square(k)) * e
vcomb = exb_part + mult * pf_part
