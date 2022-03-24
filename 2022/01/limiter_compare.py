import netCDF4
import numpy as np
import matplotlib.pyplot as plt


root = "/Users/zamperini/Documents/d3d_work/sput_testing/"
ncpaths = ["167196-sput-016.nc", "167196-sput-017a.nc", "167196-sput-017b.nc",
  "167196-sput-017c.nc", "167196-sput-017d.nc", "167196-sput-017e.nc"]
labels = ["100 m/s", "250 m/s", "500 m/s", "5 m2/s", "10 m2/s", "20 m2/s"]
pol_idx          = 21
include_neutrals = True
dist_above_lim   = 0.01
xlim             = [-2, 2]

ncs = []
for ncpath in ncpaths:
    ncs.append(netCDF4.Dataset(root + ncpath))

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

fig, ax1 = plt.subplots()
for i in range(0, len(ncs)):
    ax1.plot(xs[i][:-1], ys[i], label=labels[i])
ax1.set_xlim(xlim)
ax1.legend()
ax1.grid()
ax1.set_ylim([1e15, 5e17])
ax1.set_yscale("log")
ax1.set_xlabel("Parallel to B (m)")
ax1.set_ylabel("C Density (m-3)")
fig.tight_layout()
fig.show()
