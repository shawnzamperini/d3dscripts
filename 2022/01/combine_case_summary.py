# This script is to be used to combine two indentical 3DLIM runs differing just
# in the one uses physical sputtering and the other uses chemical.
import netCDF4
import numpy as np
import matplotlib.pyplot as plt


# Use convention that case+a = physical, case+b is chemical.
case = "022"
#ncpath_phys = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-{}a.nc".format(case)
#ncpath_chem = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-{}b.nc".format(case)
ncpath_phys = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-006.nc"
ncpath_chem = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-008.nc"
nc_phys = netCDF4.Dataset(ncpath_phys)
nc_chem = netCDF4.Dataset(ncpath_chem)

def load_output_dict(nc, dist_above_lim=0.01, include_neutrals=True,
  pol_idx=21):
    """
    Load a dictionary with a select portion of the netcdf data, as well as
    some data that has been post processed to pull out trends.

    nc (Dataset): The netCDF file.
    dist_above_lim (float): Distance above limiter to provide parallel DDLIM3
      data for.
    include_neutrals (bool): Include neutrals in DDLIM3 data.
    pol_idx (int): Poloidal index of which to pull data from.
    """

    output = {}

    # Extract some data up front.
    xs = nc.variables["XS"][:].data
    ys = nc.variables["YS"][:].data
    xwids = nc.variables["XWIDS"][:].data
    ywids = nc.variables["YWIDS"][:].data
    youts = nc.variables['YOUTS'][:].data
    xouts = nc.variables['XOUTS'][:].data
    absfac = nc.variables["ABSFAC"][:].data
    nerods3 = nc.variables["NERODS3"][:].data[4]  # 0 = deposition. 1 = primary erosion. 2 = erosion + self-sputtering. 4 = net erosion
    qtembs = nc.variables["QTEMBS"][:][0][:-2]  # Last 2 values seem to be junk.
    qtembsi = nc.variables["QTEMBSI"][:][0][:-2]
    #csintb = nc.variables["CSINTB"]
    qrnbs = nc.variables["QRNBS"][:][0][:-2]

    # Mirror the Y data around 0.
    ys = np.append(np.append(-ys[::-1], 0), ys)
    ywids = np.append(np.append(-ywids[::-1], 0), ywids)
    youts = np.append(np.append(-youts[::-1], 0), youts)

    # The longest part of the process is summing ddlim3 across the charge states.
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

    # Create Y coordinates for the DDLIM3 data above the limiter.
    y_coords = np.append(youts, youts[-1] + ywids[-1] / 2)
    above_idx = np.argmin(np.abs(xouts - dist_above_lim))
    above_ddlim3 = ddlim3[:, above_idx]

    # The net erosion at the poloidal index.
    nerods = nerods3[pol_idx] * absfac

    output["above_y"] = y_coords
    output["above_ddlim3"] = above_ddlim3
    output["nerods"] = nerods
    output["qtembs"] = qtembs
    output["qtembsi"] = qtembsi
    output["qrnbs"] = qrnbs

    return output

print("Loading physical results...")
phys = load_output_dict(nc_phys)
print("Loading chemical results...")
chem = load_output_dict(nc_chem)

# Get the maximum density at the specified distance above the limiter.
print("Density above limiter")
phys_max_nc = phys["above_ddlim3"].max()
chem_max_nc = chem["above_ddlim3"].max()
max_sum = phys_max_nc + chem_max_nc
phys_perc_str = "{:.1f}%".format(100 * phys_max_nc / max_sum)
chem_perc_str = "{:.1f}%".format(100 * chem_max_nc / max_sum)
print("{:10} {:10} {:10}".format("Physical", "Chemical", "Total"))
print("{:<10.2e} {:<10.2e} {:<10.2e}".format(phys_max_nc, chem_max_nc, max_sum))
print("{:10} {:10}".format(phys_perc_str, chem_perc_str))
print()

print("Net Erosion")
net_ero = phys["nerods"] + chem["nerods"]  # atoms/m2/s

# Where the maximum net erosion is, find how much of each it is.
idx = np.argmax(net_ero)
phys_contr = phys["nerods"][idx]
chem_contr = chem["nerods"][idx]

max_ero = net_ero.max()
tot_ero = net_ero.sum()  # atoms/m/s
phys_perc_str = "{:.1f}%".format(100 * phys_contr / net_ero.max())
chem_perc_str = "{:.1f}%".format(100 * chem_contr / net_ero.max())
print("{:10} {:10} {:10}".format("Physical", "Chemical", "Total"))
print("{:<10.2e} {:<10.2e} {:<10.2e}".format(phys["nerods"].max(), chem["nerods"].max(), net_ero.max()))
print("{:10} {:10}".format(phys_perc_str, chem_perc_str))
print()

print("Impact Energy")
eimp = 2 * phys["qtembsi"] + 3 * phys["qtembs"]
print("Max: {:.2f}".format(eimp.max()))
print("Avg: {:.2f}".format(eimp.mean()))
print()

# The deuterium flux assuming an angle of 14.
print("D Flux")
cs = 9.79e3 * np.sqrt((phys["qtembs"] + phys["qtembsi"]) / 2)
csintb = np.sin(np.radians(14))
flux = phys["qrnbs"] * cs * csintb
print("Max: {:.2e}".format(flux.max()))
print("Avg: {:.2e}".format(flux.mean()))
print()
