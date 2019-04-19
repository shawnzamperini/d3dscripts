import oedge_plots as oedge
import numpy as np
import matplotlib as mpl
import pretty_plots as pp
from scipy.optimize import curve_fit


chord_num = np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
                       14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                       28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39])

rminrsep_omp = np.array([ 0.04959334,  0.04415673,  0.03864268,  0.03157565,  0.02589525,
        0.01990257,  0.01510637,  0.01230542,  0.00897185,  0.00587476,
        0.00250091, -0.00036771, -0.00429976, -0.00667083, -0.01037882,
       -0.01331125, -0.01787746, -0.02058037, -0.02438115, -0.02737931,
       -0.03148346, -0.03533141, -0.03919977, -0.04281275, -0.0467232 ,
       -0.05064823, -0.0562748 , -0.06476167, -0.07361473, -0.08194489,
       -0.09088506, -0.10920853, -0.1416859 , -0.19046342, -0.18731829,
       -0.20910885, -0.2314993 , -0.25806067, -0.29946213, -0.32053469])

te_omp = np.array([  10.48618159,   13.30745053,   12.11449251,   14.84023123,
         15.52941079,   22.05736084,   21.10133972,   24.03670626,
         35.29908113,   43.36713562,   48.22544193,   63.17007904,
         79.35823441,   93.48981094,  130.92531738,  165.58209152,
        179.44913483,  184.20782013,  209.20890503,  268.68448181,
        249.7601532 ,  270.74040985,  289.98236084,  286.36246033,
        310.73382263,  333.06986694,  356.66596069,  388.99466858,
        433.45451355,  484.15821228,  522.49084778,  600.63817749,
        754.9612854 ,  991.8151062 , 1039.24215088, 1195.60528564,
       1333.04555664, 1525.30157471, 1763.48890381, 1869.41574707])

ne_omp = np.array([5.47127541e+18, 5.97226352e+18, 5.95883229e+18, 7.89933480e+18,
       7.60254742e+18, 8.22470498e+18, 8.38754408e+18, 9.25522274e+18,
       1.01829584e+19, 1.06047933e+19, 1.19504683e+19, 1.36506442e+19,
       1.51247046e+19, 1.65734952e+19, 1.72982687e+19, 1.86135151e+19,
       1.93561632e+19, 2.13210919e+19, 2.01411792e+19, 1.97929922e+19,
       2.14680767e+19, 2.11779494e+19, 2.17628947e+19, 2.36177393e+19,
       2.32551481e+19, 2.32474553e+19, 2.33823772e+19, 2.42153043e+19,
       2.53540558e+19, 2.45501472e+19, 2.62672233e+19, 2.77676582e+19,
       2.93389074e+19, 3.31729968e+19, 3.34515627e+19, 3.47612019e+19,
       3.65485764e+19, 3.97787412e+19, 4.24919608e+19, 4.31568698e+19])

nc_path = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/Z0-167196-004.nc'
op = oedge.OedgePlots(nc_path)

r0 = op.nc['R0'][:]
z0 = op.nc['Z0'][:]
rsep = op.rvertp[op.korpg[op.irsep-1,:op.nks[op.irsep-1]]][:,0]
zsep = op.zvertp[op.korpg[op.irsep-1,:op.nks[op.irsep-1]]][:,0]
right_half = rsep > r0
rsep = rsep[right_half]
zsep = zsep[right_half]
zidx = np.where(np.abs(zsep - z0) == np.abs(zsep-z0).min())[0][0]
rsep_omp = rsep[zidx]

r_omp = rminrsep_omp + rsep_omp
r_vals = r_omp
z_vals = np.full(len(r_vals), z0)

#op_z  = np.array([])
op_r  = np.array([])
op_te = np.array([])
op_ne = np.array([])
for i in range(0, len(chord_num)):

    # Scan through the rings.
    for ir in range(op.nrs):

        # Scan through the knots.
        for ik in range(op.nks[ir]):

            # Get the cell index of this knot on this ring.
            index = op.korpg[ir,ik] - 1

            # Only if the area of this cell is not zero append the corners.
            if op.area[ir,ik] != 0.0:
                vertices = list(zip(op.rvertp[index][0:4], op.zvertp[index][0:4]))

                # Print out a warning is the cell center is not within the vertices.
                cell = mpl.path.Path(list(vertices))
                r = op.rs[ir, ik]
                z = op.zs[ir, ik]
                if cell.contains_point([r_vals[i], z_vals[i]]):
                    print('Chord {} in cell ({}, {})'.format(i, ir, ik))
                    print('  (R, Z) Chord: ({}, {})'.format(r_vals[i], z_vals[i]))
                    print('  (R, Z) Cell:  ({}, {})'.format(r, z))
                    op_te = np.append(op_te, op.nc['KTEBS'][:][ir][ik])
                    op_ne = np.append(op_ne, op.nc['KNBS'][:][ir][ik])
                    #op_z  = np.append(op_z, z)
                    op_r  = np.append(op_r, r)

# Restrict to just outside the separatrix. Informal, but anything above Z = 0.75 is SOL here.
sol   = np.where(op_r > rsep_omp)[0]
op_r  = op_r[sol]
op_ne = op_ne[sol]
op_te = op_te[sol]

sol_ts   = np.where(r_vals > rsep_omp)[0]
r_vals   = r_vals[sol_ts]
ne_omp = ne_omp[sol_ts]
te_omp = te_omp[sol_ts]

def exp_fit(x, a, b):
    return a * np.exp(-x * b)

op_r   = op_r - rsep_omp
r_vals = r_vals - rsep_omp

popt_op, pcov_op = curve_fit(exp_fit, op_r, op_ne/1e18, maxfev=10000)
popt_ts, pcov_ts = curve_fit(exp_fit, r_vals, ne_omp/1e18, maxfev=10000)

fit_r = np.linspace(rsep_omp, r_omp.max(), 100) - rsep_omp
fit_op = exp_fit(fit_r, *popt_op)
fit_ts = exp_fit(fit_r, *popt_ts)

fig = pp.pplot(fit_r*100, fit_op*1e18, fmt='--', color=8)
fig = pp.pplot(fit_r*100, fit_ts*1e18, fmt='--', color=6, fig=fig)
fig = pp.pplot(r_vals*100, ne_omp, fmt='.', label='TS', ms=20, fig=fig)
fig = pp.pplot(op_r*100, op_ne, fmt='.', label='OEDGE', fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='ne (m-3)', color=8, ms=20)

print('OEDGE lambda ne: {:.2f} cm'.format(1/popt_op[1]*100))
print('TS lambda ne:    {:.2f} cm'.format(1/popt_ts[1]*100))
