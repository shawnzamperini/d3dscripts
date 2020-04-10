import oedge_plots
import pretty_plots as pp
import matplotlib.pyplot as plt
import numpy as np


# Load in the OedgePlots objects.
rev_path = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/Z1-167196-067.nc'
for_path = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/Z1-167196-068.nc'
rev_op = oedge_plots.OedgePlots(rev_path)
for_op = oedge_plots.OedgePlots(for_path)
rev_op.add_dat_file(rev_path.split(".nc")[0] + ".dat")
for_op.add_dat_file(for_path.split(".nc")[0] + ".dat")

def compare_along(ring=25, charge='all', fontsize=16):

    # Get the x, y Mach data.
    rev_mx, rev_my = rev_op.along_ring(ring=ring, dataname='Mach', charge=charge, plot_it=False)
    for_mx, for_my = for_op.along_ring(ring=ring, dataname='Mach', charge=charge, plot_it=False)

    # Drop last 4 (they're unecessary zeros).
    rev_mx = rev_mx[:-4]; rev_my = rev_my[:-4]
    for_mx = for_mx[:-4]; for_my = for_my[:-4]

    # Get the x, y impurity data (all charge states).
    rev_ix, rev_iy = rev_op.along_ring(ring=ring, dataname='DDLIMS', charge=charge, plot_it=False)
    for_ix, for_iy = for_op.along_ring(ring=ring, dataname='DDLIMS', charge=charge, plot_it=False)

    # Get s location of IMP, OMP, IXP and OXP for this ring.
    zs = rev_op.nc['ZS'][ring].data[:-4]  # Last four are garbage.
    zxp = float(rev_op.nc['ZXP'][ring].data)

    # IMP/OMP occur where sign of Z changes.
    zs_sign = np.sign(zs)
    signchange = ((np.roll(zs_sign, 1) - zs_sign) != 0).astype(int)
    imp_ik = np.where(signchange==1)[0][0]
    omp_ik = np.where(signchange==1)[0][1]
    imp = rev_mx[imp_ik]
    omp = rev_mx[omp_ik]

    # Same procedure, but can just do zs-zxp and find the sign changes there.
    zs_sign = np.sign(zs-zxp)
    signchange = ((np.roll(zs_sign, 1) - zs_sign) != 0).astype(int)
    ixp_ik = np.where(signchange==1)[0][0]
    oxp_ik = np.where(signchange==1)[0][1]
    ixp = rev_mx[ixp_ik]
    oxp = rev_mx[oxp_ik]

    # One way to plot it.
    """
    fig, ax1 = plt.subplots()
    ax1.plot(rev_mx, rev_my, linestyle='--', color=pp.tableau20[12])
    ax1.plot(for_mx, for_my, linestyle='-',  color=pp.tableau20[12])
    ax1.set_ylabel('Mach Number', color=pp.tableau20[12])
    ax1.set_xlabel('Along Field Line (m)')
    ax1.tick_params(axis='y', labelcolor=pp.tableau20[12])
    ax1.set_ylim([-0.3, 0.3])
    ax1.axhline(0, color='k', linestyle='--')

    ax2 = ax1.twinx()
    ax2.plot(rev_ix, rev_iy, linestyle='--', color=pp.tableau20[18], label='Reverse')
    ax2.plot(for_ix, for_iy, linestyle='-',  color=pp.tableau20[18], label='Forward')
    ax2.set_ylabel('Impurity Density (m-3)', color=pp.tableau20[18])
    ax2.tick_params(axis='y', labelcolor=pp.tableau20[18])
    ax2.set_ylim([0, 0.5e17])
    ax2.legend()

    fig.tight_layout()
    fig.show()
    """

    # A different way to plot it.
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 7))

    ax1.plot(rev_mx, rev_my, lw=3, color=pp.tableau20[8], label='Reverse')
    ax1.plot(for_mx, for_my, linestyle='-', lw=3, color=pp.tableau20[6], label='Forward')
    ax1.set_ylabel('Mach Number', fontsize=fontsize)
    #ax1.set_ylim([-0.3, 0.3])
    ax1.set_ylim([-1, 1])
    ax1.axhline(0, color='k', linestyle='-')
    ax1.axvline(imp, color='k', linestyle=(0, (5, 10)))
    ax1.axvline(omp, color='k', linestyle=(0, (5, 10)))
    ax1.axvline(ixp, color='k', linestyle=(0, (1, 10)))
    ax1.axvline(oxp, color='k', linestyle=(0, (1, 10)))
    ax1.legend(fontsize=fontsize-2, loc='upper left')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize-2)

    ax2.plot(rev_ix, rev_iy, lw=3, color=pp.tableau20[8])
    ax2.plot(for_ix, for_iy, lw=3, color=pp.tableau20[6])
    ax2.set_ylabel('Impurity Density (m-3)', fontsize=fontsize)
    ax2.set_ylim([0, 0.15e17])
    ax2.axvline(imp, color='k', linestyle=(0, (5, 10)))
    ax2.axvline(omp, color='k', linestyle=(0, (5, 10)))
    ax2.axvline(ixp, color='k', linestyle=(0, (1, 10)))
    ax2.axvline(oxp, color='k', linestyle=(0, (1, 10)))
    ax2.set_xlabel('Along Field Line (m)', fontsize=fontsize)
    ax2.tick_params(axis='both', which='major', labelsize=fontsize-2)

    fig.tight_layout()
    fig.show()


def compare_probes(match_jose=False):

    # Number convention for the fake probes is 1-5, counter clockwise starting
    # from the inner target area.
    rev_p1x, rev_p1y = rev_op.fake_probe(r_start=1.33, z_start=-1.07, r_end=1.05, z_end=-0.75, data='Mach', plot='psin', show_plot=False)
    for_p1x, for_p1y = for_op.fake_probe(r_start=1.33, z_start=-1.07, r_end=1.05, z_end=-0.75, data='Mach', plot='psin', show_plot=False)

    rev_p2x, rev_p2y = rev_op.fake_probe(r_start=1.087, z_start=0, r_end=1.02, z_end=0, data='Mach', plot='psin', show_plot=False)
    for_p2x, for_p2y = for_op.fake_probe(r_start=1.087, z_start=0, r_end=1.02, z_end=0, data='Mach', plot='psin', show_plot=False)

    rev_p3x, rev_p3y = rev_op.fake_probe(r_start=1.6, z_start=0.89, r_end=1.6, z_end=1.06, data='Mach', plot='psin', show_plot=False)
    for_p3x, for_p3y = for_op.fake_probe(r_start=1.6, z_start=0.89, r_end=1.6, z_end=1.06, data='Mach', plot='psin', show_plot=False)

    rev_p4x, rev_p4y = rev_op.fake_probe(r_start=2.25, z_start=0, r_end=2.3, z_end=0, data='Mach', plot='psin', show_plot=False)
    for_p4x, for_p4y = for_op.fake_probe(r_start=2.25, z_start=0, r_end=2.3, z_end=0, data='Mach', plot='psin', show_plot=False)

    rev_p5x, rev_p5y = rev_op.fake_probe(r_start=1.35, z_start=-1.08, r_end=1.7, z_end=-1.08, data='Mach', plot='psin', show_plot=False)
    for_p5x, for_p5y = for_op.fake_probe(r_start=1.35, z_start=-1.08, r_end=1.7, z_end=-1.08, data='Mach', plot='psin', show_plot=False)


    # Plot grid showing generally where they are. Option to match Jose's plot.
    if match_jose:
        mult = -1
        for_col = pp.tableau20[0]
        rev_col = pp.tableau20[6]
    else:
        mult = 1
        for_col = pp.tableau20[6]
        rev_col = pp.tableau20[8]

    fig, axs = plt.subplots(3, 3)
    axs[2, 0].plot(rev_p1x, mult*rev_p1y, color=rev_col)
    axs[2, 0].plot(for_p1x, mult*for_p1y, color=for_col)
    axs[1, 0].plot(rev_p2x, mult*rev_p2y, color=rev_col)
    axs[1, 0].plot(for_p2x, mult*for_p2y, color=for_col)
    axs[0, 1].plot(rev_p3x, mult*rev_p3y, color=rev_col)
    axs[0, 1].plot(for_p3x, mult*for_p3y, color=for_col)
    axs[1, 2].plot(rev_p4x, mult*rev_p4y, color=rev_col)
    axs[1, 2].plot(for_p4x, mult*for_p4y, color=for_col)
    axs[2, 2].plot(rev_p5x, mult*rev_p5y, color=rev_col)
    axs[2, 2].plot(for_p5x, mult*for_p5y, color=for_col)

    # Add horizontal line at 0, plot limits.
    for ax in axs.flatten():
        #ax.set_ylim([-0.5, 0.5])
        ax.set_ylim([-1.0, 1.0])
        ax.set_xlim([1.0, 1.10])
        ax.axhline(0.0, color='k', linestyle='--')

    # For the right plots invert the x axis.
    axs[2, 0].set_xlim([1.10, 1.0])
    axs[1, 0].set_xlim([1.10, 1.0])
    axs[0, 0].set_xlim([1.10, 1.0])

    fig.tight_layout()
    fig.show()
