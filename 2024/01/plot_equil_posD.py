#.............................................#
#.
#.Plot geometry from gkyl grid file.
#.Manaure Francisquez.
#.
#.............................................#
import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from mpl_toolkits.mplot3d import Axes3D
import sys
#sys.path.insert(0, '/Users/manaure/Documents/codebits/pets/gkeyll/postProcessingScripts/')
import pgkylUtil_adios2 as pgu
from freeqdsk import geqdsk
from scipy.interpolate import CubicSpline


dataDir = ['/scratch/gpfs/tbernard/CEDA/TCV_shots/posD/']
outDir  = dataDir

#.Root name of files to process.
simName = ['tcv-axisym-posD-3x2v']

Raxis   = 0.86 #0.8727315068 #.Major radius at the magnetic axis (m).

#.Indicate what you want to plot below.
plot3D            = False
plotPoloidalPlane = False

outFigureFile      = True     #.Output a figure file?.
figureFileFormat   = '.png'    #.Can be .png, .pdf, .ps, .eps, .svg.
outDir             = './'      #.Output directory where you wish to save figures.

polyOrder = 1
basis     = 'ms'
basisName = "Ser"


# ..................... end of user inputs (MAYBE) ......................... #

#.Some fontsizes used in plots.
xyLabelFontSize       = 17
titleFontSize         = 17
colorBarLabelFontSize = 17
tickFontSize          = 14
legendFontSize        = 14


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Read in eqdsk files
# Positive triangularity

# # Retrive gkyl grid values
# fI = 0
# fileRoot = dataDir[fI] + simName[fI]
# gridNodes = pgu.getRawData(fileRoot + '_grid.bp')
# Rp, Zp, zp = list(), list(), list()
#
# simDim = np.size(np.shape(gridNodes))-1
# slices = [np.s_[0:1]]*(simDim+1)
# slices[simDim] = 0
# for d in range(simDim):
#   slices[d] = np.s_[0:]
# numCells = np.shape(gridNodes[tuple(slices)])
#
# nodeList = np.reshape(gridNodes,[np.product(numCells),3])
#
# #.Compute field line following coordinate z to color the nodes.
# #.Below we refer to Cartesian coordinates (X,Y,Z), cylindrical
# #.coordinates (R,phi,Z) and toroidal coordinates (r,theta,phi).
# #.The field aligned coordinate is z.
# X, Y, Z = nodeList[:,0], nodeList[:,1], nodeList[:,2]
# phi   = np.arctan2(Y,X)
# R     = X/np.cos(phi)
# theta = np.arctan2(Z,(R-Raxis))
# r     = (R-Raxis)/np.cos(theta) #  = R*cos(phi)
# eps   = r/Raxis
#
# Rp.append(R)
# Zp.append(Z)

# retrieve exp
with open("default/g065125.01298","r") as f:   # edit accordingly to gfile
  data = geqdsk.read(f)

psi_RZ = data["psi"]
qpsi = data["qpsi"]
psi_1d = np.linspace(data["simagx"], data["sibdry"], data["nx"])

R = np.linspace(data["rleft"], data["rleft"]+data["rdim"], data["nx"])
Z = np.linspace(data["zmid"] - data["zdim"]/2, data["zmid"] + data["zdim"]/2, data["ny"])
RR, ZZ = np.meshgrid(R, Z, indexing='ij')

print('Raxis =', data["rmagx"])
print('Zaxis =', data["zmagx"])

Rlim = np.amax(data["rlim"]) - np.amin(data["rlim"])
Zlim = np.amax(data["zlim"]) - np.amin(data["zlim"])
levels = np.linspace(psi_RZ.min(), psi_RZ.max(), 20)

# find indices in psi to map to R
Zmag_diff = np.abs(Z - data["zmagx"])
zmag_idx = np.where(Zmag_diff==np.amin(Zmag_diff))
print('Z index at Baxis = ', zmag_idx)
Zidx = int(input('Z index? '))
psiR_zmag = psi_RZ[:,Zidx]
Rmag_idx = np.where(psiR_zmag==np.amin(psiR_zmag))
print('R index at Baxis = ', Rmag_idx)
Ridx = int(input('R index? '))
psiR_zmag = psiR_zmag[Ridx:]

# map psi value to R
RofPsi = CubicSpline(psiR_zmag, R[62:])
Rforq = RofPsi(psi_1d)

# Plot experimental equilibria 
fig = plt.figure(figsize=(Rlim*5+2,Zlim*5))
ax = fig.add_subplot(111)
ax.plot(data["rlim"],data["zlim"])
ax.plot(data["rbdry"],data["zbdry"])
ax.set_xlabel('R (m)')
ax.set_ylabel('Z (m)')
psi_con = ax.contourf(RR, ZZ, psi_RZ,levels=levels)
fig.colorbar(psi_con)
# for fI in range(int(0.66*len(Rp[0]))):
#   ax.scatter(Rp[0][fI], Zp[0][fI],color='black',s=1)
# for fI in range(int(0.66*len(Rp[0])),len(Rp[0])):
#   ax.scatter(Rp[0][fI], Zp[0][fI],color='red',s=1)
ax.set_aspect("equal")
plt.tight_layout()
if outFigureFile:
  plt.savefig(dataDir[0] + simName[0] + '_gridPointsCartesian' + figureFileFormat)
plt.show()
plt.close()

Rforq = Rforq[len(Rforq)//2:]
qpsi = qpsi[len(qpsi)//2:]
qPoly = np.polyfit(Rforq, qpsi,3)
print(qPoly)
Rpoly = np.linspace(Rforq[0], Rforq[-1], 100)
qProf = qPoly[3] + qPoly[2]*Rpoly + qPoly[1]*Rpoly**2 + qPoly[0]*Rpoly**3

plt.plot(Rforq, qpsi, label='exp')
plt.plot(Rpoly, qProf, label='poly')
plt.legend()
plt.xlabel('R (m)')
plt.ylabel('q')
plt.show()
plt.close()

# pressure
# pr = data["pres"]
# pr = pr[len(pr)//2:]
# plt.plot(Rforq,pr/2e19/1.602e-19)
# plt.show()
# plt.close()

print('LCFS =', Rforq[-1])

