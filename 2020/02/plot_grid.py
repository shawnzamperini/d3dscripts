import netCDF4
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


path = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/gridtest_d3d_carre.nc'
nc = netCDF4.Dataset(path)

# Create a mesh of of the corners of the each cell/polygon in the grid.
#mesh  = np.array([])
mesh = []
num_cells = 0
rs     = nc['RS'][:]
zs     = nc['ZS'][:]
nrs    = nc['NRS'][:]
nks    = nc['NKS'][:]
korpg  = nc['KORPG'][:]
rvertp = nc['RVERTP'][:]
zvertp = nc['ZVERTP'][:]
rvesm  = nc['RVESM'][:]
zvesm  = nc['ZVESM'][:]
# Scan through the rings.
for ir in range(nrs):

    # Scan through the knots.
    for ik in range(nks[ir]):

        # Get the cell index of this knot on this ring.
        index = korpg[ir,ik] - 1

        # Only if the area of this cell is not zero append the corners.
        #if self.area[ir,ik] != 0.0:
        vertices = list(zip(rvertp[index][0:4], zvertp[index][0:4]))
        #mesh = np.append(mesh, vertices)
        mesh.append(vertices)
        num_cells = num_cells + 1

        # Print out a warning is the cell center is not within the vertices.
        cell = mpl.path.Path(list(vertices))
        r = rs[ir, ik]
        z = zs[ir, ik]
        if not cell.contains_point([r, z]):
            print("Error: Cell center not within vertices.")
            print("  (ir, ik)    = ({}, {})".format(ir, ik))
            print("  Vertices    = {}".format(vertices))
            print("  Cell center = ({}, {})".format(r, z))


fig = plt.figure()
ax = fig.add_subplot(111)
coll = mpl.collections.PolyCollection(mesh)
ax.add_collection(coll)
ax.plot(rvesm, zvesm, color='k', linewidth=1)
ax.axis('equal')
ax.set_xlim([0.9,2.5])
