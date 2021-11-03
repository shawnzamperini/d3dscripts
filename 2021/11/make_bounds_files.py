# Quick helper script to create a set of .bounds files for 3DLIM.
import sys
import numpy as np
from LimWallToolkit import LimWallToolkit


# OPTION 2 FROM lwt_control_file.py.
lwt = LimWallToolkit()
for tor in np.arange(0, 360, 30, dtype=int):
    print(tor)

    # These files are from running MAFOT in the +1 and -1 directions. Which
    # corresponds to file1 and to file2 depends on your setup. One way to make it
    # work is file1 can contains all the connection lengths in the positive
    # 3DLIM direction (i.e. that which replaces surface 1, L19 in the 3DLIM input
    # file), and file2 be the negative values (surface 2, L21).
    mafot_file1 = "/Users/zamperini/Documents/d3d_work/184527/mafot/lam_hires_tor{}_conn+1.dat".format(tor)
    mafot_file2 = "/Users/zamperini/Documents/d3d_work/184527/mafot/lam_hires_tor{}_conn-1.dat".format(tor)

    toroidal_angle = tor  # MiMES = 240.

    # The R bins for 3DLIM. Can have them already set or set here and copy/paste
    # into 3DLIM input file.
    lim_rbins = np.arange(-0.1, 0.02, 0.0025)

    # The number of Pbins MUST match 2*MAXNPS+1 in 3DLIM, otherwise it will not
    # work. Probably equal to 41. Also I think it needs to be symmetric...
    lim_pbins = np.linspace(-0.7, 0.7, 41)

    # In machine coordinates, where do you want the origin of 3DLIM to be?
    # MiMES-like origin.
    r_origin2 = 2.295
    z_origin2 = -0.188
    # DiMES-like origin.
    #r_origin2 = 1.485
    #z_origin2 = -1.15

    # Which machine coordinate does the 3DLIM "R" coordinate go along?
    # MiMES = "R", DiMES = "Z" and None will set it along the actual R-Rsep
    # direction (i.e. the plasma radial direction).
    along_coord = "R"

    # Must have .bound extension, i.e. 184527_conns.bound or ramp.bound are okay.
    output_bounds_file = "/Users/zamperini/Documents/lim_runs/184527_tor{}.bound".format(tor)
    gfile_pickle_path2 = "/Users/zamperini/Documents/d3d_work/184527/184527_3500.pickle"
    wall_path2 = "/Users/zamperini/Documents/d3d_work/184527/mafot/mafot_3D_wall.dat"

    dbg = lwt.bounds_file_from_mafot(tor_angle=toroidal_angle,
      mafot_file1=mafot_file1, mafot_file2=mafot_file2, lim_rbins=lim_rbins,
      lim_pbins=lim_pbins, r_origin=r_origin2,
      z_origin=z_origin2, output_file=output_bounds_file,
      gfile_pickle_path=gfile_pickle_path2, along_coord=along_coord,
      wall_path=wall_path2, show_plot=False)
