import numpy as np
import matplotlib.pyplot as plt
import netCDF4


ncpath  = "/Users/zamperini/Documents/lim_runs/colprobe-a8-actual-001.nc"
nc = netCDF4.Dataset(ncpath)

def get_shit(nc, vary=True):

    # Pull out relevant arrays.
    ys = nc.variables["YS"][:].data
    xs = nc.variables["XS"][:].data
    ps = nc.variables["PS"][:].data
    ywids = nc.variables["YWIDS"][:].data
    xwids = nc.variables["XWIDS"][:].data
    pwids = nc.variables["PWIDS"][:].data
    if vary:
        vp3d = nc.variables["velplasma_4d_2"][:].data
        #vp3d = nc.variables["velplasma_4d_2"][:].data
    else:
        vp3d = nc.variables["velplasma"][1,:,:].data
        #vp3d = nc.variables["velplasma"][1,:,:].data
    ddlim3 = nc.variables["DDLIM3"][:].sum(axis=1)

    # Go through appropriately named files and get averages of DDLIM3.
    #count = 1
    #for letter in ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"]:
    #    try:
    #        tmp_nc = netCDF4.Dataset(ncpath.split(".nc")[0] + letter + ".nc")
    #        print("Loaded {} netcdf file".format(letter))
    #        tmp_ddlim3 = tmp_nc.variables["DDLIM3"][:].sum(axis=1)
    #        ddlim3 += tmp_ddlim3
    #        count += 1
    #    except:
    #        pass

    # Average ddlim3.
    #ddlim3 = ddlim3 / count

    # Calculate bin centers.
    rad_locs = xs - xwids / 2
    pol_locs = ps - pwids / 2

    # Special treatment for the Y coordinate since it needs to be mirrored and
    # a zero added.
    tmp = ys - ywids / 2
    par_locs = np.append(np.append(-tmp[::-1], 0), tmp)

    # Identify where unecessary zeros are. DO this way to avoid accidentally
    # removing the zero that has been added in the previous step.
    y_keep_start = np.nonzero(par_locs)[0].min()
    y_keep_end   = np.nonzero(par_locs)[0].max()
    x_keep_start = np.nonzero(rad_locs)[0].min()
    x_keep_end   = np.nonzero(rad_locs)[0].max()
    #p_keep_start = np.nonzero(pol_locs)[0].min()
    #p_keep_end   = np.nonzero(pol_locs)[0].max()

    # Index all variables accordingly.
    rad_locs = rad_locs[x_keep_start:x_keep_end]
    #pol_locs = pol_locs[p_keep_start:p_keep_end]
    par_locs = par_locs[y_keep_start:y_keep_end]
    #te3d = te3d[y_keep_start:y_keep_end, x_keep_start:x_keep_end, p_keep_start:p_keep_end]
    #te3d = te3d[y_keep_start:y_keep_end, x_keep_start:x_keep_end]
    ddlim3 = ddlim3[:, y_keep_start:y_keep_end, x_keep_start:x_keep_end]
    vp3d = vp3d[y_keep_start:y_keep_end, x_keep_start:x_keep_end]

    # Mask the zeros for a better plot.
    #te3d = np.ma.masked_where(te3d == 0, te3d)
    ddlim3 = np.ma.masked_where(ddlim3 == 0, ddlim3)

    # Not to be efficient, so just lazy copy/paste.
    dep_arr = np.array(nc.variables['NERODS3'][0] * -1)
    cp_ps     = np.array(nc.variables['PS'][:].data)
    cp_pwids  = np.array(nc.variables['PWIDS'][:].data)
    cp_pol_locs = ps - pwids/2.0
    cp_rad_locs = np.array(nc.variables['ODOUTS'][:].data)
    cline = np.abs(pol_locs).min()

    # Index the deposition array at the centerline for plotting.
    itf_x = cp_rad_locs[np.where(cp_rad_locs > 0.0)[0]]
    itf_y = dep_arr[np.where(cp_pol_locs == cline)[0], np.where(cp_rad_locs > 0.0)[0]]
    otf_x = cp_rad_locs[np.where(cp_rad_locs < 0.0)[0]] * -1
    otf_y = dep_arr[np.where(cp_pol_locs == cline)[0], np.where(cp_rad_locs < 0.0)[0]]

    return {"par_locs":par_locs, "rad_locs":rad_locs, "pol_loc":pol_locs,
      "ddlim3":ddlim3, "vp3d":vp3d, "itf_x":itf_x, "itf_y":itf_y, "otf_x":otf_x,
      "otf_y":otf_y, "dep_arr":dep_arr}

cp = get_shit(nc)

dl3_sum = cp["ddlim3"].sum(axis=0).T
dep = cp["dep_arr"]

fig, axs = plt.subplots(1, 3, figsize=(8, 5), sharex=True, sharey=True)
axs[0].set_xlim([-6, 6])
axs[0].set_facecolor("grey")
axs[1].set_facecolor("grey")
axs[2].set_facecolor("grey")

c1 = axs[0].pcolormesh(cp["par_locs"], cp["rad_locs"], dl3_sum, shading="auto", vmin=0, vmax=0.5, cmap="magma")

fig.colorbar(c1, ax=axs[0])

fig.tight_layout()
fig.show()
