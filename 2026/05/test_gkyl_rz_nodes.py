# This script is adapted from Antoine's script

import postgkyl as pg
import numpy as np
from scipy.interpolate import RegularGridInterpolator

def custom_meshgrid(x,y,z=0):
    # custom meshgrid function to have natural orientation (x,y,z)
    if np.isscalar(z):
        if len(x.shape) > 1 and len(y.shape) > 1:
            return [x,y]
        Y,X = np.meshgrid(y,x)
        return [X,Y]
    else:
        if len(x.shape) > 1 and len(y.shape) > 1 and len(z.shape) > 1:
            return [x,y,z]
        Y,X,Z = np.meshgrid(y,x,z)
        return [X,Y,Z]

field_frame = Frame(self.sim, fieldname=fieldName, tf=timeFrame, load=True)
self.gridsN = field_frame.xNodal # nodal grids
self.ndim = len(self.gridsN) # Dimensionality


#.Centered mesh creation
ndim = 3
meshC = [[] for i in range(ndim)] 
for i in range(ndim):
	nNodes  = len(self.gridsN[i])
	meshC[i] = np.zeros(nNodes-1)
	meshC[i] = np.multiply(0.5,self.gridsN[i][0:nNodes-1]+self.gridsN[i][1:nNodes])
self.dimsC = [np.size(meshC[i]) for i in range(self.ndim)]
self.meshC = meshC

xxI, zzI = math_tools.custom_meshgrid(self.meshC[0],self.zgridI)
self.dimsI = np.shape(xxI) # interpolation plane dimensions (R,Z)

# Path to the nodes file
root = "/mnt/d/gkeyll_results/west-sol-62104"
nodes_path = root + "/gk_west_lsn_sol_3x2v_p1-nodes.gkyl"


# Load nodes. Copying procedure from here: https://github.com/Antoinehoff/
#   personal_gkyl_scripts/blob/add_neut_species/pygkyl/pygkyl/projections/
#   poloidalprojection.py#L295
nodalData = pg.GData(nodes_path)
nodalVals = nodalData.get_values()
alpha_idx = 0
R = nodalVals[:, alpha_idx, :, 0]
Z = nodalVals[:, alpha_idx, :, 1]
Phi = nodalVals[:, alpha_idx, :, 2] 

nodalGridTemp = nodalData.get_grid()   # contains one more element than number of nodes.
nodalGrid = []
for d in range(0,len(nodalGridTemp)):
	nodalGrid.append(np.linspace(nodalGridTemp[d][0], nodalGridTemp[d][-1], len(nodalGridTemp[d])-1))

RInterpolator = RegularGridInterpolator((nodalGrid[0], nodalGrid[2]), R)
ZInterpolator = RegularGridInterpolator((nodalGrid[0], nodalGrid[2]), Z)
PhiInterpolator = RegularGridInterpolator((nodalGrid[0], nodalGrid[2]), Phi)

Rint = RInterpolator((xxI, zzI))
Zint = ZInterpolator((xxI, zzI))
Phiint = PhiInterpolator((xxI, zzI))

      if self.sim.geom_param.geom_type == 'efit':
        self.alpha_rz_phi0 = -self.gridsN[1][alpha_idx] - Phiint # Overwrite the results of compute_alpha
        phiTor = np.pi
        for k in range(self.kyDimsC[1]): # Overwrite the results of compute_xyz2RZ
          for iz in range(self.nzI):
            shift = -2*np.pi*self.alpha_rz_phi0[:,iz]/self.LyC + phiTor
            self.xyz2RZ[:,+k,iz]  = np.exp(1j*k*shift)
            #.Negative ky's.
            self.xyz2RZ[:,-k,iz] = np.conj(self.xyz2RZ[:,+k,iz])

    self.RIntN, self.ZIntN = np.zeros((self.dimsI[0]+1,self.dimsI[1]+1)), np.zeros((self.dimsI[0]+1,self.dimsI[1]+1))
    for j in range(self.dimsI[1]):
        for i in range(self.dimsI[0]):
            self.RIntN[i,j] = Rint[i,j]-0.5*(Rint[1,j]-Rint[0,j])
        self.RIntN[self.dimsI[0],j] = Rint[-1,j]+0.5*(Rint[-1,j]-Rint[-2,j])
        self.RIntN[:,self.dimsI[1]] = self.RIntN[:,-2]

    for i in range(self.dimsI[0]):
        for j in range(self.dimsI[1]):
            self.ZIntN[i,j] = Zint[i,j]-0.5*(Zint[i,1]-Zint[i,0])
        self.ZIntN[i,self.dimsI[1]] = Zint[i,-1]+0.5*(Zint[i,-1]-Zint[i,-2])
        self.ZIntN[self.dimsI[0],:] = self.ZIntN[-2,:]
