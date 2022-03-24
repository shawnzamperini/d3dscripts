import netCDF4
import matplotlib.pyplot as plt
import numpy as np


root = "/Users/zamperini/Documents/d3d_work/sput_testing/"
c_dens = [0, 1, 2, 3]
ti_mult = [2, 3, 4, 5]
pol_idx = 21

ncs = {}
for c in c_dens:
    for m in ti_mult:
        path = "{}167196-5cm-lim-{}-{}.nc".format(root, m, c)
        ncs["{}-{}".format(m, c)] = netCDF4.Dataset(path)

# A reference case.
nc_3_25 = netCDF4.Dataset("{}/167196-5cm-lim-3-25.nc".format(root))

def get_nerods_data(nc):

    output = {}

    # Extract some data up front.
    xs = nc.variables["XS"][:].data
    ys = nc.variables["YS"][:].data
    xwids = nc.variables["XWIDS"][:].data
    ywids = nc.variables["YWIDS"][:].data
    youts = nc.variables['YOUTS'][:].data
    xouts = nc.variables['XOUTS'][:].data
    absfac = nc.variables["ABSFAC"][:].data
    qtembs = nc.variables["QTEMBS"][:][0][:-2]  # Last 2 values seem to be junk.
    qtembsi = nc.variables["QTEMBSI"][:][0][:-2]
    cthetb = float(nc.variables["CTHETB"][:].data)
    qrnbs = nc.variables["QRNBS"][:][0][:-2]
    odouts = nc.variables["ODOUTS"][:]
    cl = float(nc['CL'][:].data)
    ca  = float(nc['CA'][:].data)
    caw = float(nc['CAW'][:].data)

    # The limiter locations.
    qxs = nc.variables["QXS"][:].data
    qedges1 = nc.variables["QEDGES"][:].data[0]
    qedges2 = nc.variables["QEDGES"][:].data[1]

    # Erosion data. 0 = deposition. 1 = primary erosion.
    # 2 = erosion + self-sputtering. 4 = net erosion
    nerods3 = nc.variables["NERODS3"][:].data[2][pol_idx] * absfac

    # Only need the values up to len(qedges). Trim zeros. Want one zero at the
    # end for fill between.
    nonzero = np.nonzero(qedges1)
    qxs = qxs[:len(qedges1)][nonzero]
    qedges1 = qedges1[nonzero]
    qedges2 = qedges2[nonzero]
    qxs = np.append(qxs, 0)
    qedges1 = np.append(qedges1, 0)
    qedges2 = np.append(qedges2, 0)

    # Calculate the surface flux as calculated in the code. Need to mirror it
    # so it aligns with odouts and makes sense.
    cs = 9.79e3 * np.sqrt((qtembs + qtembsi) / 2)
    csintb = np.sin(np.radians(cthetb))
    flux = qrnbs * cs * csintb
    half = int(len(flux)/2)
    flux = np.append(flux[half:], flux[half:][::-1])

    #output["X"] = X
    #output["Y"] = Y
    #output["ddlim3"] = ddlim3
    #output["above_y"] = y_coords
    #output["above_ddlim3"] = above_ddlim3
    output["odouts"] = odouts
    output["nerods3"] = nerods3
    output["qtembs"] = qtembs
    output["qtembsi"] = qtembsi
    output["qrnbs"] = qrnbs
    output["qedges1"] = qedges1
    output["qedges2"] = qedges2
    output["qxs"] = qxs
    output["ca"] = ca
    output["caw"] = caw
    output["cl"] = cl
    output["flux"] = flux

    return output

# Pull out data for each case.
outputs = {}
for key, nc in ncs.items():
    outputs[key] = get_nerods_data(nc)
output_3_25 = get_nerods_data(nc_3_25)

# Assemble into a 2D array for plotting.
Y, X = np.meshgrid(ti_mult, c_dens)
Z = np.zeros(X.shape)
for ci in range(0, len(c_dens)):
    for mi in range(0, len(ti_mult)):

        # Straight up maximum gross erosion.
        key = "{}-{}".format(ti_mult[mi], c_dens[ci])
        ero = outputs[key]["nerods3"].max()

        Z[mi, ci] = ero


# 2D heatmap of the erosion data.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

cont = ax1.contourf(X, Y, Z)
ax1.scatter(X.flatten(), Y.flatten(), c="k", s=50)
cbar = fig.colorbar(cont, ax=ax1)
ax1.set_ylabel("Ti Multiplier")
ax1.set_xlabel("Background C (%ne)")

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, len(ti_mult)))
for ci in range(0, len(c_dens)):
#    color = colors[ci]
#for mi in range(0, len(ti_mult)):
    color = colors[ci]
    label = "{}%".format(c_dens[ci])
    key = "{}-{}".format(2, c_dens[ci])
    x = outputs[key]["odouts"]
    y = outputs[key]["nerods3"]
    ax2.plot(x, y, color=color, label=label)
ax2.plot(output_3_25["odouts"], output_3_25["nerods3"], label="25%", color="r")
ax2.legend()

fig.tight_layout()
fig.show()
