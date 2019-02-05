import numpy as np
import matplotlib.pyplot as plt
import netCDF4


def itf_over_otf(netcdf):
    """

    """
    #The deposition array.
    dep_arr = np.array(netcdf.variables['NERODS3'][0] * -1)

    # Location of each P bin, and its width. Currently they all have the same width,
    # but it may end up such that there are custom widths so we leave it like this.
    ps     = np.array(netcdf.variables['PS'][:].data)
    pwids  = np.array(netcdf.variables['PWIDS'][:].data)

    # Array of poloidal locations (i.e. the center of each P bin).
    pol_locs = ps - pwids/2.0

    # Drop last row since it's garbage.
    dep_arr = dep_arr[:-1, :]
    pol_locs = pol_locs[:-1]

    # Distance cell centers along surface (i.e. the radial locations).
    rad_locs = np.array(netcdf.variables['ODOUTS'][:].data)

    # Get the centerline index (or closest to it).
    cline = np.abs(pol_locs).min()

    # Index the deposition array at the centerline for plotting.
    itf_x = rad_locs[np.where(rad_locs > 0.0)[0]]
    itf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
    otf_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
    otf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]

    itf_div_otf = itf_y / otf_y

    return itf_x, itf_div_otf


fig = plt.figure()
ax = fig.add_subplot(111)

for test in ['24', '25', '26', '27', '28', '29', '30', '32', '33', '34']:
    filename = '/home/shawn/Documents/d3d_work/3DLIM Runs/colprobe-z0-0' + test + '.nc'
    netcdf = netCDF4.Dataset(filename)
    x, y = itf_over_otf(netcdf)
    ax.plot(x, y, '.', label=test)

ax.set_xlabel('Distance along probe (m)')
ax.set_ylabel('ITF/OTF')
ax.legend()
fig.show()
