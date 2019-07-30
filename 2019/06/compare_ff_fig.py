import oedge_plots
import numpy as np


# Some constants used.
zero_sub = 1e-20
vz_mult = 0.9

print("Loadin' OEDGE data...")

# Load OedgePlots object.
ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/Z0-167196-001.nc'
op = oedge_plots.OedgePlots(ncpath)

print("Loadin' forces...")

# Load forces we want to compare.
ff5  = op.calculate_forces('FF', charge=5, vz_mult=vz_mult).data
ff10 = op.calculate_forces('FF', charge=10, vz_mult=vz_mult).data
ff15 = op.calculate_forces('FF', charge=15, vz_mult=vz_mult).data
ff20 = op.calculate_forces('FF', charge=20, vz_mult=vz_mult).data
ff30 = op.calculate_forces('FF', charge=30, vz_mult=vz_mult).data
ff40 = op.calculate_forces('FF', charge=40, vz_mult=vz_mult).data
fig5   = op.calculate_forces('FIG', charge=5)
fig10  = op.calculate_forces('FIG', charge=10)
fig15  = op.calculate_forces('FIG', charge=15)
fig20  = op.calculate_forces('FIG', charge=20)
fig30  = op.calculate_forces('FIG', charge=30)
fig40  = op.calculate_forces('FIG', charge=40)

print("Messin' with the data...")

# We only want to look at the region where the net forces are upstream, which on
# the right half is > 0, and left half is < 0. Thus we need to know which half we
# are on so we can filter out the right data. The middle should be where
# KSS > KSMAXS / 2.0.
kss    = op.nc['KSS'][:]
ksmaxs = op.nc['KSMAXS'][:]
right_half = kss > ksmaxs[:, np.newaxis] / 2.0
raw_data = right_half

# Let's put this 2D array into the same format as the 1D array, using the same
# method in read_data_2d.
data = np.zeros(op.num_cells)

count = 0
for ir in range(op.nrs):
    for ik in range(op.nks[ir]):
        if op.area[ir, ik] != 0.0:
            data[count] = raw_data[ir][ik]
            count = count + 1

# Great. Now let's say if the corresponding force is right_half == 0 (so in the
# left half), then multiply by -1. This will just make sure a positive value
# always means an upstream force, and can easily grab only those later.
for force in [ff5, ff10, ff15, ff20, ff30, ff40, fig5, fig10, fig15, fig20, fig30, fig40]:
    for i in range(0, len(data)):
        if data[i] == 0:
            force[i] *= -1

        # So now let's just discard all the downstream forces (the positive ones)
        # by replacing them with zeros.
        if force[i] > 0:
            force[i] = 0.0

        # Zeros goof up the colorbar, so just make them a real small number as a
        # quick fix.
        if force[i] == 0.0:
            force[i] = zero_sub

# Plot a couple to make sure they look right.
#op.plot_contour_polygon('FF10', own_data=ff10, normtype='symlog', vmin=-1e-15, vmax=1e-15)
#op.plot_contour_polygon('FIG',  own_data=fig,  normtype='symlog', vmin=-1e-15, vmax=1e-15)

# Now create the ratios. If < 0, FF dominates. If > 0, FIG dominates.
ratio5  = np.zeros(len(fig5))
ratio10 = np.zeros(len(fig5))
ratio15 = np.zeros(len(fig5))
ratio20 = np.zeros(len(fig5))
ratio30 = np.zeros(len(fig5))
ratio40 = np.zeros(len(fig5))
for i in range(0, len(ratio10)):
    if ff5[i]  != zero_sub and fig5[i] != zero_sub:
        ratio5[i]  = fig5[i]/ff5[i]  - 1.0
    if ff10[i] != zero_sub and fig10[i] != zero_sub:
        ratio10[i] = fig10[i]/ff10[i] - 1.0
    if ff15[i] != zero_sub and fig15[i] != zero_sub:
        ratio15[i] = fig15[i]/ff15[i] - 1.0
    if ff20[i] != zero_sub and fig20[i] != zero_sub:
        ratio20[i] = fig20[i]/ff20[i] - 1.0
    if ff30[i] != zero_sub and fig30[i] != zero_sub:
        ratio30[i] = fig30[i]/ff30[i] - 1.0
    if ff40[i] != zero_sub and fig40[i] != zero_sub:
        ratio40[i] = fig40[i]/ff40[i] - 1.0

print("Makin' some plots...")

# Plot the ratio.
v = 2
op.plot_contour_polygon('<--- FF FiG --->',  own_data=ratio5,  normtype='symlin', vmin=-v, vmax=v)
#op.plot_contour_polygon('<--- FF10+ FiG --->', own_data=ratio10, normtype='symlin', vmin=-v, vmax=v)
#op.plot_contour_polygon('<--- FF15+ FiG --->', own_data=ratio15, normtype='symlin', vmin=-v, vmax=v)
#op.plot_contour_polygon('<--- FF20+ FiG --->', own_data=ratio20, normtype='symlin', vmin=-v, vmax=v)
#op.plot_contour_polygon('<--- FF30+ FiG --->', own_data=ratio30, normtype='symlin', vmin=-v, vmax=v)
#op.plot_contour_polygon('<--- FF40+ FiG --->', own_data=ratio40, normtype='symlin', vmin=-v, vmax=v)
