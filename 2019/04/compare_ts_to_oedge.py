import oedge_plots as oedge
import numpy as np
import matplotlib as mpl
import pretty_plots as pp
from scipy.optimize import curve_fit


chord_num = np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
                       14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                       28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39])
avg_dens = np.array([4.4517439e+18, 5.3568102e+18, 5.5170739e+18, 7.0700120e+18,
                     7.3630836e+18, 7.7147855e+18, 8.3650097e+18, 8.9819677e+18,
                     8.9858418e+18, 9.5908970e+18, 1.1757347e+19, 1.2785762e+19,
                     1.5182987e+19, 1.6761432e+19, 1.8515730e+19, 1.8030272e+19,
                     1.9498563e+19, 2.1502757e+19, 1.9889090e+19, 2.0171179e+19,
                     2.1396144e+19, 2.1486639e+19, 2.2232077e+19, 2.2914293e+19,
                     2.2764064e+19, 2.3631872e+19, 2.3717484e+19, 2.4286620e+19,
                     2.5334492e+19, 2.4899681e+19, 2.6708299e+19, 2.7809165e+19,
                     2.9627291e+19, 3.3348751e+19, 3.3569906e+19, 3.5530600e+19,
                     3.6975846e+19, 3.9987325e+19, 4.2539542e+19, 4.3528214e+19])
avg_temp = np.array([   8.812114,   11.245648,   11.508032,   13.754086,   16.479189,
                       19.329979,   22.973291,   26.130997,   33.522972,   39.07105 ,
                       42.090603,   54.697388,   77.85025 ,  100.078026,  126.938255,
                      161.58322 ,  173.65654 ,  179.61266 ,  209.85075 ,  243.1496  ,
                      232.14796 ,  257.38324 ,  268.5137  ,  288.8754  ,  320.15622 ,
                      322.3992  ,  357.149   ,  392.24805 ,  437.68353 ,  479.59534 ,
                      521.12384 ,  594.9005  ,  742.7704  ,  995.41254 , 1018.54645 ,
                     1146.3402  , 1295.6997  , 1525.4973  , 1765.0964  , 1865.0049])
r_vals = np.array([1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999,
                   1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999,
                   1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999,
                   1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999,
                   1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999,
                   1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999, 1.9399999,
                   1.9399999, 1.9399999, 1.9399999, 1.9399999])
z_vals = np.array([ 0.8253    ,  0.8138    ,  0.80230004,  0.78779995,  0.77629995,
                    0.7643    ,  0.7548    ,  0.7493    ,  0.7428    ,  0.7368    ,
                    0.73029995,  0.7248    ,  0.7173    ,  0.71279997,  0.7058    ,
                    0.7003    ,  0.6918    ,  0.6868    ,  0.67980003,  0.6743    ,
                    0.66679996,  0.6598    ,  0.65279996,  0.64629996,  0.6393    ,
                    0.63229996,  0.62229997,  0.6073    ,  0.5918    ,  0.5773    ,
                    0.5618    ,  0.53029996,  0.4748    ,  0.3908    ,  0.3963    ,
                    0.35779998,  0.3168    ,  0.26479998,  0.16479999,  0.0268])

nc_path = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/Z0-167196-004.nc'
op = oedge.OedgePlots(nc_path)

op_z  = np.array([])
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
                    op_z  = np.append(op_z, z)
                    op_r  = np.append(op_r, r)

# Restrict to just outside the separatrix. Informal, but anything above Z = 0.75 is SOL here.
sol   = np.where(op_z > 0.75)[0]
op_z  = op_z[sol]
op_ne = op_ne[sol]
op_te = op_te[sol]

sol_ts   = np.where(z_vals > 0.75)[0]
z_vals   = z_vals[sol_ts]
avg_dens = avg_dens[sol_ts]
avg_temp = avg_temp[sol_ts]

def exp_fit(x, a, b):
    return a * np.exp(-x * b)

popt_op, pcov_op = curve_fit(exp_fit, op_z, op_ne/1e18, maxfev=10000)
popt_ts, pcov_ts = curve_fit(exp_fit, z_vals, avg_dens/1e18, maxfev=10000)

fit_z = np.linspace(0.75, 0.85, 100)
fit_op = exp_fit(fit_z, *popt_op)
fit_ts = exp_fit(fit_z, *popt_ts)

fig = pp.pplot(fit_z, fit_op*1e18, fmt='--', color=8)
fig = pp.pplot(fit_z, fit_ts*1e18, fmt='--', color=6, fig=fig)
fig = pp.pplot(z_vals, avg_dens, fmt='.', label='TS', ms=20, fig=fig)
fig = pp.pplot(op_z, op_ne, fmt='.', label='OEDGE', fig=fig, xlabel='Z (m)', ylabel='ne (m-3)', color=8, ms=20)

print('OEDGE lambda ne: {:.2f} cm'.format(1/popt_op[1]*100))
print('TS lambda ne:    {:.2f} cm'.format(1/popt_ts[1]*100))
