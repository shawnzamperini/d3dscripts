import oedge_plots
import numpy as np
import matplotlib.pyplot as plt
import pretty_plots as pp


plt.rcParams['font.family'] = 'DejaVu Sans'

# The input files we want.
nc3 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/archive/Z0-167196-003.nc'
nc64 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/Z0-167196-064.nc'

# Constants
charge = 10
vz_mult = 0.99
ring = 20

op3  = oedge_plots.OedgePlots(nc3)
op64 = oedge_plots.OedgePlots(nc64)

ffx3,   ffy3   = op3.along_ring(ring,  'ff', charge=charge, vz_mult=vz_mult, plot_it=False)
figx3,  figy3  = op3.along_ring(ring,  'fig', charge=charge, plot_it=False)
fegx3,  fegy3  = op3.along_ring(ring,  'feg', charge=charge, plot_it=False)
fex3,   fey3   = op3.along_ring(ring,  'fe', charge=charge, plot_it=False)
ffx64,  ffy64  = op64.along_ring(ring, 'ff', charge=charge, vz_mult=vz_mult, plot_it=False)
figx64, figy64 = op64.along_ring(ring, 'fig', charge=charge, plot_it=False)
fegx64, fegy64 = op64.along_ring(ring, 'feg', charge=charge, plot_it=False)
fex64,  fey64  = op64.along_ring(ring, 'fe', charge=charge, plot_it=False)

"""
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ffx3[:-4], ffy3[:-4], '-', label='FF', color=pp.tableau20[8])
ax.plot(figx3[:-4], figy3[:-4], '-', label='FiG', color=pp.tableau20[6])
#ax.plot(fegx3, np.abs(fegy3), 'b-')
#ax.plot(fex3, np.abs(fey3), 'k-')
ax.set_xlabel('S (m)')
ax.set_ylabel('Force (N)')
ax.set_yscale('symlog', linthreshy=1e-21)
fig.tight_layout()

fnet = ffy3[:-4]+figy3[:-4]
fig = pp.pplot(ffx3[:-4],  ffy3[:-4], fmt='-', label='FF', color=20, linestyle='dotted')
fig = pp.pplot(ffx3[:-4],  figy3[:-4], fmt='-', label='FiG', color=20, fig=fig, linestyle='dashed')
fig = pp.pplot(figx3[:-4][net_idx], fnet[net_idx], fmt='-', label='Fnet', fig=fig, xlabel='s (m)',
               ylabel='Force (N)', logy=True, linthreshy=1e-20, color=20, weight='bold')

# Want a matching font for the slides.
plt.rcParams['font.family'] = 'DejaVu Sans'


# Want grid lines.
yticks = np.array([1e-16, 1e-17, 1e-18, 1e-19, 1e-20])
fig.axes[0].set_yticks(np.append(yticks, -yticks), minor=True)
fig.axes[0].grid(True, which='minor')
fig.axes[0].grid(True, which='major')

fig.show()
"""

# Want to get rid of the craziness when we're close to zero.
#net_idx = np.append(np.arange(0, 78), np.arange(133, 256))
net_idx = np.append(np.arange(0, 88), np.arange(123, 256))

def do_plot(ffx, ffy, figx, figy):
    fnet = ffy[:-4]+figy[:-4]
    fig = pp.pplot(ffx[:-4],  ffy[:-4], fmt='-', label='FF', color=20, linestyle='dotted')
    fig = pp.pplot(ffx[:-4],  figy[:-4], fmt='-', label='FiG', color=20, fig=fig, linestyle='dashed')
    fig = pp.pplot(figx[:-4][net_idx], fnet[net_idx], fmt='-', label='Fnet', fig=fig, xlabel='s (m)',
                   ylabel='Force (N)', logy=True, linthreshy=1e-20, color=20, weight='bold')

    # Want a matching font for the slides.
    plt.rcParams['font.family'] = 'DejaVu Sans'


    # Want grid lines.
    yticks = np.array([1e-16, 1e-17, 1e-18, 1e-19, 1e-20])
    fig.axes[0].set_yticks(np.append(yticks, -yticks), minor=True)
    fig.axes[0].grid(True, which='minor')
    fig.axes[0].grid(True, which='major')

    fig.show()

do_plot(ffx3, ffy3, figx3, figy3)
do_plot(ffx64, ffy64, figx64, figy64)
