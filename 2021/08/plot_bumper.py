import numpy as np
import matplotlib.pyplot as plt
import netCDF4


# Load in netcdf file.
ncpath  = "/Users/zamperini/Documents/lim_runs/varying-wall-016.nc"
ncpath_no  = "/Users/zamperini/Documents/lim_runs/no-varying-wall-015.nc"
nc = netCDF4.Dataset(ncpath)
nc_no = netCDF4.Dataset(ncpath_no)


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

    return {"par_locs":par_locs, "rad_locs":rad_locs, "pol_loc":pol_locs,
      "ddlim3":ddlim3, "vp3d":vp3d}

vary = get_shit(nc, True)
no_vary = get_shit(nc_no, False)

dl3_sum_vary = vary["ddlim3"].sum(axis=0).T
dl3_sum_no = no_vary["ddlim3"].sum(axis=0).T
vp3d_mid_vary = vary["vp3d"][:,:,25].T
vp3d_mid_no = no_vary["vp3d"].T

# Three poloidal plots.
ip1 = 7
ip2 = 15
ip3 = 35

#Z1 = te3d[:,:,ip1].T
#Z2 = te3d[:,:,ip2].T
#Z3 = te3d[:,:,ip3].T

# Sum over a poloidal range for better stats.
Z1 = dl3_sum_vary; vmin1=0; vmax1=0.005; cmap="magma"
Z2 = dl3_sum_no; tit1 = "Varying bound"; tit2 = "Normal bound"; tit3 = "Vary / Normal"
Z3 = dl3_sum_vary/dl3_sum_no; vmin2 = -2; vmax2=2

#Z1 = vp3d_mid_vary; vmin1=-10000; vmax1=10000; cmap="coolwarm"
#Z2 = vp3d_mid_no; tit1 = "Varying bound"; tit2 = "Normal bound"; tit3 = "Vary / Normal"
#Z3 = vp3d_mid_vary/vp3d_mid_no; vmin2 = -5; vmax2=5

#Z1 = vary["ddlim3"][5,:,:].T; vmin1=0; vmax1=0.00015; cmap="magma"
#Z2 = vary["ddlim3"][23,:,:].T; tit1 = "Ramp #1"; tit2 = "Ramp #2"; tit3 = "Normal"
#Z3 = no_vary["ddlim3"][23,:,:].T; vmin2 = -0; vmax2=0.00015

# Some plots.
fig, axs = plt.subplots(1, 3, figsize=(8, 5), sharex=True, sharey=True)
axs[0].set_xlim([-5, 5])
axs[0].set_facecolor("grey")
axs[1].set_facecolor("grey")
axs[2].set_facecolor("grey")

#cmap = "magma"
c1 = axs[0].pcolormesh(vary["par_locs"], vary["rad_locs"], Z1, shading="auto", vmin=vmin1, vmax=vmax1, cmap=cmap)
c2 = axs[1].pcolormesh(no_vary["par_locs"], no_vary["rad_locs"], Z2, shading="auto", vmin=vmin1, vmax=vmax1, cmap=cmap)
c3 = axs[2].pcolormesh(vary["par_locs"], vary["rad_locs"], Z3, shading="auto", cmap=cmap, vmin=vmin2, vmax=vmax2)
#axs[2].pcolormesh(par_locs, rad_locs, Z3, shading="auto", cmap="magma")

fig.colorbar(c1, ax=axs[0])
fig.colorbar(c2, ax=axs[1])
fig.colorbar(c3, ax=axs[2])

#axs[0].set_title("P = {:.2f}".format(pol_locs[ip1]))
#axs[1].set_title("P = {:.2f}".format(pol_locs[ip2]))
#axs[2].set_title("P = {:.2f}".format(pol_locs[ip3]))
axs[0].set_title(tit1)
axs[1].set_title(tit2)
axs[2].set_title(tit3)

fig.tight_layout()
fig.show()
