from read_gfile import read_gfile
import matplotlib.pyplot as plt
import numpy as np
import flan_plots
from skimage import measure


# First load gfile for WEST #62104 and pull out the R, Z grid and the 2D psi
# array defined on it
path = "west_LSN.geqdsk"
gfile = read_gfile(path)
rgrid = gfile["rgrid"]
zgrid = gfile["zgrid"]
psirz = gfile["psirz"]

# Then load Flan simulation
fp = flan_plots.FlanPlots("/pscratch/sd/z/zamp/flandir/west_nearsol_point/ot_point_source.nc")

# Pull out the grid coordinates
x = fp.nc["geometry"]["x"][:]  # psi
y = fp.nc["geometry"]["y"][:]  # alpha
z = fp.nc["geometry"]["z"][:]  # chi or theta

# Pull out R, Z contours of each flux tube, defined by psi value
for psi in x:
	contours = measure.find_contours(psirz, psi)
