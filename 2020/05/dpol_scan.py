import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import lim_plots
from scipy.interpolate import interp1d

plt.rcParams['font.family'] = 'DejaVu Sans'
# These are the "Tableau 20" colors as RGB. The last one is just black. I added it.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
             (0, 0, 0)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

# LAMS file to compare against.
#cu04_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
#            'Polodial_Scans/New Map Script Results/CU04_Map_Analysis.xlsx'
#cd04_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
#            'Polodial_Scans/New Map Script Results/CD04_Map_Analysis.xlsx'
cu04_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/CD07_Map_Analysis.xlsx'
cd04_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/CU07_Map_Analysis.xlsx'

# The ncfiles of our Dpol scan.
ncpath0_20 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-049e.nc'  # Dpol = 0.2
ncpath0_50 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-049f.nc'  # Dpol = 0.5
ncpath2_00 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-049g.nc'  # Dpol = 2.0
ncpath0_05 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-049h.nc'  # Dpol = 0.05
ncpath0_02 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-049i.nc'  # Dpol = 0.02
ncpath0_005 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-049j.nc'  # Dpol = 0.005

# Get average poloidal profile from LAMS for first 5 cm.
cu04 = pd.read_excel(cu04_file).pivot_table(columns='z Location [mm]', index='Axial Location [mm]')
cd04 = pd.read_excel(cd04_file).pivot_table(columns='z Location [mm]', index='Axial Location [mm]')
rad_cutoff = 50
zloc_otf = cu04.iloc[cu04.index<=rad_cutoff].mean()['Total W'].index
zloc_itf = cd04.iloc[cd04.index<=rad_cutoff].mean()['Total W'].index
otf_totw = cu04.iloc[cu04.index<=rad_cutoff].mean()['Total W'].values
itf_totw = cd04.iloc[cd04.index<=rad_cutoff].mean()['Total W'].values
otf_totw_err = cu04.iloc[cu04.index<=rad_cutoff].std()['Total W'].values
itf_totw_err = cd04.iloc[cd04.index<=rad_cutoff].std()['Total W'].values

# Get the average poloidal profile from LAMS data where at each "slice" of a
# bathub curve it is normalized to the max. Then each "slice" is averaged.
bathtubs = len(cu04.iloc[cu04.index<=rad_cutoff].index)
all_tubs_otf = np.zeros((bathtubs, len(cu04.iloc[0]['Total W'].index)))
all_tubs_itf = np.zeros((bathtubs, len(cd04.iloc[0]['Total W'].index)))
for i in range(0, bathtubs):

    # Normalized bathtub at each radial slice.
    norm_tub_otf = cu04.iloc[i]['Total W'] / cu04.iloc[i]['Total W'].max()
    norm_tub_itf = cd04.iloc[i]['Total W'] / cd04.iloc[i]['Total W'].max()
    all_tubs_otf[i] = norm_tub_otf
    all_tubs_itf[i] = norm_tub_itf

# Then the average normalized bathtub.
avg_tub_otf = all_tubs_otf.mean(axis=0)
std_tub_otf = all_tubs_otf.std(axis=0)
avg_tub_itf = all_tubs_itf.mean(axis=0)
std_tub_itf = all_tubs_itf.std(axis=0)

# Get average poloidal profiles likewise from 3DLIM.
def get_avg_pol(ncpath):
    """
    Function to get the average poloidal data from a 3DLIM run, or a combination
    of them.
    """
    rad_cutoff = 0.05
    #rad_cutoff = 0.10

    # Load in the deposition array.
    nc = netCDF4.Dataset(ncpath)
    print(ncpath)
    try:
        dep_arr = np.array(nc.variables['NERODS3'][0] * -1)
    except:
        dep_arr = np.zeros((6, 2*nc.variables['MAXNPS'][:]+1, nc.variables['MAXOS'][:]))
        print("  No NERODS3.")

    # Add on contributions from repeat runs.
    for i in range(1, 11):
        try:
            #print('Loading {}...'.format(i))
            ncpath_add = ncpath.split('.nc')[0] + str(i) + '.nc'
            nc = netCDF4.Dataset(ncpath_add)
            print("Found additional run: {}".format(ncpath_add))
            try:
                dep_arr = dep_arr + np.array(nc.variables['NERODS3'][0] * -1)
            except KeyError:
                print("  No NERODS3.")
        except:
            pass

    # Code copied from lim_plots..
    dep_arr = np.array(nc.variables['NERODS3'][0] * -1)
    ps     = np.array(nc.variables['PS'][:].data)
    pwids  = np.array(nc.variables['PWIDS'][:].data)
    pol_locs = ps - pwids/2.0
    dep_arr = dep_arr[:-1, :]
    pol_locs = pol_locs[:-1]
    rad_locs = np.array(nc.variables['ODOUTS'][:].data)
    idx = np.where(np.abs(rad_locs)<rad_cutoff)[0]
    rad_locs = rad_locs[idx]
    dep_arr = dep_arr[:, idx]
    idx = np.where(rad_locs > 0.0)[0]
    X_itf, Y_itf = np.meshgrid(rad_locs[idx], pol_locs)
    #Z_itf = dep_arr[:, idx]

    # Normalize the data on a per row (i.e. per radial slice) basis.
    Z_itf = dep_arr[:, idx][:, ::-1] / dep_arr[:, idx][:, ::-1].max(axis=0)

    idx = np.where(rad_locs < 0.0)[0]
    X_otf, Y_otf = np.meshgrid(np.abs(rad_locs[idx][::-1]), pol_locs)
    #Z_otf = dep_arr[:, idx][:, ::-1]

    # Normalize the data on a per row (i.e. per radial slice) basis.
    max_vals = dep_arr[:, idx][:, ::-1].max(axis=0)
    z_vals = dep_arr[:, idx][:, ::-1]
    bad_idxs = np.where(max_vals == 0)[0]
    if len(bad_idxs) > 0:
        print("Warning: Statistics are bad. Consider more runs.")
        max_vals = np.delete(max_vals, bad_idxs)
        z_vals = np.delete(z_vals, bad_idxs, 1)
    #Z_otf = dep_arr[:, idx][:, ::-1] / dep_arr[:, idx][:, ::-1].max(axis=0)
    Z_otf = z_vals / max_vals
    #print('OTF Max: {}'.format(dep_arr[:, idx][:, ::-1].max(axis=0)))

    # Average along the radial direction.
    avg_pol_itf = np.mean(Z_itf, 1)
    avg_pol_otf = np.mean(Z_otf, 1)

    # Get the centerline index (or closest to it).
    cline = np.abs(pol_locs).min()
    cline_idx = np.where(pol_locs == cline)[0][0]

    # Get average peaking factor for each side.
    peak1 = avg_pol_itf[:cline_idx].max() / avg_pol_itf[cline_idx]
    peak2 = avg_pol_itf[cline_idx:].max() / avg_pol_itf[cline_idx]
    itf_peak = (peak1 + peak2) / 2.0
    peak1 = avg_pol_otf[:cline_idx].max() / avg_pol_otf[cline_idx]
    peak2 = avg_pol_otf[cline_idx:].max() / avg_pol_otf[cline_idx]
    otf_peak = (peak1 + peak2) / 2.0

    message = "OTF/ITF Peaking Ratio: {:.2f}".format(otf_peak/itf_peak)
    print(message)

    # Shift so it starts at zero (C probe size).
    pol_locs = pol_locs + 0.0025 / 2
    return {'pol_locs':pol_locs, 'avg_pol_itf':avg_pol_itf, 'avg_pol_otf':avg_pol_otf}

d0_005 = get_avg_pol(ncpath0_005)
d0_05 = get_avg_pol(ncpath0_05)
d0_02 = get_avg_pol(ncpath0_02)
d0_20 = get_avg_pol(ncpath0_20)
d0_50 = get_avg_pol(ncpath0_50)
d2_00 = get_avg_pol(ncpath2_00)

# Restrict the data to just the x range where probe data exists, or else
# the normalization to compare real vs 3dlim doesn't make sense.
def overlap(d, side):

    x = d['pol_locs'] * 1000
    if side == 'otf':
        f_y = interp1d(x, d['avg_pol_otf'])
        x_int = np.linspace(zloc_otf.min(), zloc_otf.max(), 25)
    elif side == 'itf':
        f_y = interp1d(x, d['avg_pol_itf'])
        x_int = np.linspace(zloc_itf.min(), zloc_itf.max(), 25)
    y_int = f_y(x_int) / f_y(x_int).max()

    return x_int, y_int

x0_005, y0_005 = overlap(d0_005, 'otf')
x0_05, y0_05 = overlap(d0_05, 'otf')
x0_02, y0_02 = overlap(d0_02, 'otf')
x0_20, y0_20 = overlap(d0_20, 'otf')
x0_50, y0_50 = overlap(d0_50, 'otf')
x2_00, y2_00 = overlap(d2_00, 'otf')
x0_005_itf, y0_005_itf = overlap(d0_005, 'itf')
x0_05_itf, y0_05_itf   = overlap(d0_05, 'itf')
x0_02_itf, y0_02_itf   = overlap(d0_02, 'itf')
x0_20_itf, y0_20_itf   = overlap(d0_20, 'itf')
x0_50_itf, y0_50_itf   = overlap(d0_50, 'itf')
x2_00_itf, y2_00_itf   = overlap(d2_00, 'itf')

# Plot it.
fig, ax = plt.subplots()
ax.fill_between((zloc_otf-1.25)/10, avg_tub_otf-std_tub_otf, avg_tub_otf+std_tub_otf, alpha=0.3, color=tableau20[18])
ax.plot((zloc_otf-1.25)/10, avg_tub_otf, marker='.', linestyle='-', label='LAMS', color=tableau20[18])
#ax.plot(x0_02, y0_02, color=tableau20[12], alpha=1.0, label='3DLIM')
#ax.plot(x0_05, y0_05, color=tableau20[12], alpha=0.85)
#ax.plot(x0_20, y0_20, color=tableau20[12], alpha=0.70)
#ax.plot(x0_50, y0_50, color=tableau20[12], alpha=0.55)
#ax.plot(x2_00, y2_00, color=tableau20[12], alpha=0.4)
#ax.plot(x0_02, y0_02, color=tableau20[12], alpha=1.0, label='3DLIM')
ax.plot((x0_05-1.25)/10, y0_05, color=tableau20[12], alpha=1.0, label="3DLIM")
ax.plot((x0_20-1.25)/10, y0_20, color=tableau20[12], alpha=0.66)
ax.plot((x0_50-1.25)/10, y0_50, color=tableau20[12], alpha=0.33)
#ax.plot(x2_00, y2_00, color=tableau20[12], alpha=0.4)
ax.set_xlabel('Poloidal (cm)', fontsize=16)
ax.set_ylabel('W Deposition (normalized)', fontsize=16)
ax.set_xlim([(-0.5-1.25)/10, (3-1.25)/10])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelsize=12)
ax.legend(fontsize=16, loc='lower left')


"""
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(zloc_itf, avg_tub_itf, marker='.', linestyle='-', label='LAMS', color=tableau20[18])
ax2.fill_between(zloc_otf, avg_tub_otf-std_tub_otf, avg_tub_otf+std_tub_otf, alpha=0.3, color=tableau20[18])
ax2.plot(zloc_otf, avg_tub_otf, marker='.', linestyle='-', label='LAMS', color=tableau20[18])
ax1.plot(x0_02_itf, y0_02_itf, color=tableau20[12], alpha=1.0, label='3DLIM')
ax1.plot(x0_05_itf, y0_05_itf, color=tableau20[12], alpha=0.85)
ax1.plot(x0_20_itf, y0_20_itf, color=tableau20[12], alpha=0.70)
ax1.plot(x0_50_itf, y0_50_itf, color=tableau20[12], alpha=0.55)
ax1.plot(x2_00_itf, y2_00_itf, color=tableau20[12], alpha=0.4)
ax2.plot(x0_02, y0_02, color=tableau20[12], alpha=1.0, label='3DLIM')
ax2.plot(x0_05, y0_05, color=tableau20[12], alpha=0.85)
ax2.plot(x0_20, y0_20, color=tableau20[12], alpha=0.70)
ax2.plot(x0_50, y0_50, color=tableau20[12], alpha=0.55)
ax2.plot(x2_00, y2_00, color=tableau20[12], alpha=0.4)
ax1.set_xlabel('Poloidal (mm)', fontsize=16)
ax1.set_ylabel('Deposition (normalized)', fontsize=16)
ax2.set_xlabel('Poloidal (mm)', fontsize=16)
ax2.set_ylabel('Deposition (normalized)', fontsize=16)
ax2.set_xlim([-0.5, 3])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.tick_params(labelsize=12)
ax2.legend(fontsize=12, loc='lower left')
"""

fig.tight_layout()
fig.show()
