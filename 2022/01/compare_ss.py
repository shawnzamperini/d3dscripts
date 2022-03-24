# Script to compare C densities above a semi-circle limiter for each combination
# of physical and chemical sputtering both w and w/o self-sputtering.
import netCDF4
import numpy as np
import matplotlib.pyplot as plt


ncpath_phys_ss = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-019a.nc"
ncpath_chem_ss = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-019b.nc"
ncpath_phys_no = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-019c.nc"
ncpath_chem_no = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-019d.nc"

nc_phys_ss = netCDF4.Dataset(ncpath_phys_ss)
nc_chem_ss = netCDF4.Dataset(ncpath_chem_ss)
nc_phys_no = netCDF4.Dataset(ncpath_phys_no)
nc_chem_no = netCDF4.Dataset(ncpath_chem_no)

ncs = [nc_phys_ss, nc_chem_ss, nc_phys_no, nc_chem_no]
labels = ["Phys - SS", "Chem - SS", "Phys", "Chem"]

pol_idx          = 21
include_neutrals = True
dist_above_lim   = 0.01
xlim             = [-2, 2]

def get_par_imp(nc, dist_above_lim):

    # Extract some data up front.
    xs = nc.variables["XS"][:].data
    ys = nc.variables["YS"][:].data
    xwids = nc.variables["XWIDS"][:].data
    ywids = nc.variables["YWIDS"][:].data
    youts = nc.variables['YOUTS'][:].data
    xouts = nc.variables['XOUTS'][:].data
    absfac = nc.variables["ABSFAC"][:].data

    # Mirror the Y data around 0.
    ys = np.append(np.append(-ys[::-1], 0), ys)
    ywids = np.append(np.append(-ywids[::-1], 0), ywids)
    youts = np.append(np.append(-youts[::-1], 0), youts)

    # Extract DDLIM3, sum across all charge states, ignore neutrals if asked to.
    if include_neutrals:
        ddlim3 = nc.variables["DDLIM3"][:].data[pol_idx, 1:, :, :].sum(axis=0)
    else:
        ddlim3 = nc.variables["DDLIM3"][:].data[pol_idx, 2:, :, :].sum(axis=0)
    ddlim3 = np.ma.masked_where(ddlim3<=0, ddlim3)

    # Scale to m-3.
    ddlim3 = ddlim3 * absfac

    # We then need to trim all the trailing zeros on the arrays.
    xkeep_min = np.nonzero(xs)[0].min()
    xkeep_max = np.nonzero(xs)[0].max()
    ykeep_min = np.nonzero(ys)[0].min()
    ykeep_max = np.nonzero(ys)[0].max()
    ys     = ys[ykeep_min:ykeep_max]
    ywids  = ywids[ykeep_min:ykeep_max]
    youts  = youts[ykeep_min:ykeep_max]
    ddlim3 = ddlim3[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

    y_coords = np.append(youts, youts[-1] + ywids[-1] / 2)
    above_idx = np.argmin(np.abs(xouts - dist_above_lim))
    print("above_idx = {}, X = {}".format(above_idx, xouts[above_idx]))
    above_ddlim3 = ddlim3[:, above_idx]

    return y_coords, above_ddlim3

xs = []; ys = []
for i in range(0, len(ncs)):
    print("{}/{}".format(i+1, len(ncs)))
    x, y = get_par_imp(ncs[i], dist_above_lim)
    xs.append(x)
    ys.append(y)

fig, ax1 = plt.subplots(figsize=(5,4))

ax1.plot(xs[0][:-1], ys[0], label=labels[0], linestyle="--", color="tab:red")
ax1.plot(xs[1][:-1], ys[1], label=labels[1], linestyle="-.", color="tab:red")
ax1.plot(xs[0][:-1], ys[0]+ys[1], label="Total - SS", linestyle="-", color="tab:red")
ax1.plot(xs[2][:-1], ys[2], label=labels[2], linestyle="--", color="tab:purple")
ax1.plot(xs[3][:-1], ys[3], label=labels[3], linestyle="-.", color="tab:purple")
ax1.plot(xs[2][:-1], ys[2]+ys[3], label="Total", linestyle="-", color="tab:purple")

ax1.set_xlim([-3, 3])
ax1.set_ylim([1e15, 1e19])
ax1.legend()
ax1.set_yscale("log")
ax1.set_xlabel("Distance from limiter center (m)")
ax1.set_ylabel("C Density (m-3)")
ax1.grid()

fig.tight_layout()
fig.show()
