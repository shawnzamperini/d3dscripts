############################################
#                                          #
#  This script reads in data from NSTX SOL #
#  simulation over a user defined time     #
#  interval. Pgkyl interpolates the data   #
#  and returns 3D numpy arrays.            #
#                                          #
############################################

import postgkyl as pg

import pylab
import math
import numpy as np
import csv
import os

from matplotlib import rcParams

print('Starting data analysis...')
fstart = int(input('fstart? '))
fend = int(input('fend? '))
step = 1
Nt = fend-fstart
Nt //= step
dt = 1e-6

fp = "nstx-neut"

# physical constants
mp = 1.672623e-27
AMU = 2.014 # Deuterium ions
mi = mp*AMU
me = 9.10938188e-31
eV = 1.602e-19

c_s = np.sqrt(40*eV/mi)

# For equilibrium profiles
d = pg.GData("%s_electron_GkM0_%d.bp" %(fp,fstart))
dg = pg.GInterpModal(d, 1, 'ms')
X, nElc = dg.interpolate(0)
nElc = nElc[:,:,:,0]

zi = len(X[2])//2 # to take a slice at midplane

Nx = len(X[0])-1
Nz = len(X[2])-1
Ny = len(X[1])-1

# calculate grids
x = X[0]
dx = np.diff(x)[0]
x = x + dx/2
x = x[:-1]

y = X[1]
dy = np.diff(y)[0]
y = y + dy/2
y = y[:-1]

z = X[2]
dz = np.diff(z)[0]
z = z + dz/2
z = z[:-1]

# Create arrays for equilibrium profiles that are averaged in y, z and time
nElc_tot = []
tElc_tot = []
tIon_tot = []
phi_tot = []

for t in range(fstart,fend, step):
    if t % 100 == 0:
        print('frame = ', t)
    if (t == fend - 1):
        print("frame = ", t)
            
    # electron density
    d = pg.GData("%s_electron_GkM0_%d.bp" %(fp,t))
    dg = pg.GInterpModal(d, 1, 'ms')
    X, nElc = dg.interpolate(0)
    nElc = nElc[:,:,:,0]    
    nElc_tot.append(nElc)
    
    # electron temperature
    d = pg.GData("%s_electron_GkTemp_%d.bp" %(fp,t))
    dg = pg.GInterpModal(d, 1, 'ms')
    X, tElc = dg.interpolate(0)
    tElc = tElc[:,:,:,0]/eV
    tElc_tot.append(tElc)
    
    # ion temperature
    d = pg.GData("%s_ion_GkTemp_%d.bp" %(fp,t))
    dg = pg.GInterpModal(d, 1, 'ms')
    X, tIon = dg.interpolate(0)
    tIon = tIon[:,:,:,0]/eV
    tIon_tot.append(tIon)
    
    # Phi eq profile calcuation
    d = pg.GData("%s_phi_%d.bp" %(fp,t))
    dg = pg.GInterpModal(d, 1, 'ms')
    X, phi = dg.interpolate(0)
    phiZmin = phi[:,:,0,0]
    phi = phi[:,:,:,0]
    phi_tot.append(phi)
    
nElc_tot = np.array(nElc_tot)
tElc_tot = np.array(tElc_tot)
tIon_tot = np.array(tIon_tot)
phi_tot = np.array(phi_tot)

print("Numpy array size (Nt, Nx, Ny, Nz) = ", np.shape(nElc_tot))
