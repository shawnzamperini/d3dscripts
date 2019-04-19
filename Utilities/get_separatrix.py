"""
returns the contents of a separatrix file as a line collection.
separatrix file is copied manually from the appropriate parts of the 
extended grid file, and formatted to one line per cell.
alternatively, use get_sep_nc to parse separatrix directly from nc object.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc

def get_sep(sepfile=None):

    if not sepfile:
        sepfile='C:\Users\Jacob\Documents\separatrix_grid_167196_3000_ext.dat' # D3D 167196
    comment='*'
	
	# extract r2,z2,r4,z4 (inner sides of cells)
    coords=np.genfromtxt(sepfile,comments=comment,usecols=(5,6,9,10),unpack=True)
	
    lines=[]
    for i in range(len(coords[0])):
        lines.append([(coords[0][i],coords[1][i]),(coords[2][i],coords[3][i])])
    return lines

def get_sep_nc(ncobj=None):
    # ncobj is a loaded oedgeNC/oedgeNCimp/oedgeNCmap object
    if not ncobj:
        print 'get_separatrix.get_sep_nc: needs a loaded oedgeNC object'
        return
    
    rsep = ncobj.rvertp[ncobj.korpg[ncobj.irsep-1,:ncobj.nks[ncobj.irsep-1]]][:,0]
    zsep = ncobj.zvertp[ncobj.korpg[ncobj.irsep-1,:ncobj.nks[ncobj.irsep-1]]][:,0]
    nsep = len(rsep)
    lines=[]
    for i in range(nsep-2): # dont connect final point to first
        lines.append([(rsep[i],zsep[i]),(rsep[i+1],zsep[i+1])])
    return lines

if __name__ == '__main__':
	a=get_sep()
	#print a
	lc =mc.LineCollection(a,linewidths=1)
	ax = plt.axes()
	ax.add_collection(lc)
	ax.set_aspect('equal','datalim','C')
	ax.set_xlim([-0.2,1.8])
	ax.set_ylim([-2,2])
	plt.show()
