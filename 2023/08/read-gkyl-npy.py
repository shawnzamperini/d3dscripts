##############################################
#                                            #
#  Script reads numpy files output from      # 
#  Gkeyll 'clopen' DIII-D simulation with    #
#  postivie triangularity (delta = +0.4).    #
#                                            #
##############################################

import numpy as np

# These files have shape (Nt, Nx, Ny, Nz) = (653, 192, 192, 32)
nElc = np.load('d3d_elc_density.npy')
tElc = np.load('d3d_elc_temp.npy')
tIon = np.load('d3d_ion_temp.npy')
phi = np.load('d3d_phi.npy')

# Time array
tGrid = np.linspace(0.0, 652e-6, 653)

# Grid in field-line following coordinate system
xGrid = np.linspace(0., 0.15, 192)
yGrid = np.linspace(-7.074822e-02, 7.074822e-02, 192)
zGrid = np.linspace(-np.pi, np.pi, 32)

