# This script takes the connection length data from each toroidal angle, as
# calulated by MAFOT, and converts it into a format for input to 3DLIM via
# the .bounds extension. This mostly just means playing the 3D trig game.
import pandas as pd
import numpy  as np
import matplotlib.pyplot as plt
from mafot_read_3d_wall import read_wall
from shapely.geometry.polygon import Polygon
from matplotlib.colors import LogNorm, Normalize


# Constants used in this script.
root = "/Users/zamperini/Documents/d3d_work/184527/mafot/"
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin", "psimax",
  "psiav", "pitch angle", "yaw angle", "theta", "psi"]
max_tor = 200  # Eventually set and leave at 359.
lim_rbins = np.arange(-0.1, 0.02, 0.0025)  # Lower range of each bin.
pmin = -0.5; pmax = 0.5; pstep = 0.03
#lim_pbins = np.arange(pmin, pmax+pstep, pstep)
lim_pbins = np.linspace(pmin, pmax, 41)
r_origin = 2.29  # The machine coordinate of the R = 0 coordinate in 3DLIM.
z_origin = -0.188  # The machine coordinate of the P = 0 coordinate in 3DLIM.
tor_origin = 0  # The machine toroidal angle that the origin is at. MiMES = 240.
output_fname = "184527.bound"

# Load in all the MAFOT data for each toroidal angle.
Ls = {}; pitches = {}
first_run = True
for tor in range(0, max_tor+1):

    # Load data into DataFrame.
    fname = root + "lam_torang{}.dat".format(tor)
    df = pd.read_csv(fname, skiprows=52, names=columns, delimiter="\t")

    # Reshape into 2D arrays. The R, Z locations don't change among the
    # runs, so just load and shape them once.
    if first_run:
        r = df["R (m)"].unique()
        z = df["Z (m)"].unique()
        R, Z = np.meshgrid(r, z)
        first_run = False
    l = df["Lconn (km)"].values * 1000  # km to m
    p = df["pitch angle"].values
    L = l.reshape(len(r), len(z))
    pitch = p.reshape(len(r), len(z))
    Ls[tor] = L
    pitches[tor] = pitch

# Get the corresponding R values from the 3DLIM R bins.
rs = [r_origin - rb for rb in lim_rbins]

# The corresponding Z values from the 3DLIM P (for perpendicular) bins are a
# bit more works since the P direction is not exactly in the Z direction and is
# angled off by the pitch angle of the field line. This means the P direction
# does not stay at a single toroidal angle, and instead slightly jumps across
# different angles.
# Pitch angle at the origin (minus bc of mafot). We are approximating that the
# pitch angle does not change along the field line, otherwise a rectangular
# grid like 3DLIM does not work. Square into a circle kinda thing.
dist = np.sqrt(np.square(R - r_origin) + np.square(Z - z_origin))
nearest_idx = np.where(dist == np.nanmin(dist))
origin_pitch = pitches[tor_origin][nearest_idx]
zs = [float(z_origin + pb * np.cos(origin_pitch)) for pb in lim_pbins]

# Do a linear fit so in later scripts we can convert between lim_pbins and Z.
# Return coefficients here are highest power first.
pfit = np.polyfit(lim_pbins, zs, 1)

# Create empty bounds array and fill it with the closest connection length.
# The approximation above, where it is assumed the pitch angle does not change
# across toroidal angles is still made here.
bounds = np.zeros((len(rs), len(zs)))
L = Ls[tor_origin]
for ir in range(0, len(rs)):
    for iz in range(0, len(zs)):
        dist = np.sqrt(np.square(R - rs[ir]) + np.square(Z - zs[iz]))
        nearest_idx = np.where(dist == np.nanmin(dist))
        bounds[ir][iz] = L[nearest_idx]

# Print to a text file ready for input to 3DLIM.
with open(output_fname, "w") as f:
    print("Writing to {}".format(output_fname))
    f.write("Absorbing boundary data. Rows are for each radial bin, columns are for each poloidal bin. Created via bounds_file_from_mafot.py.\n")
    f.write("Dimensions: {} {}\n".format(len(rs), len(zs)))
    for ir in range(0, len(rs)):
        f.write("{:5.2f}".format(bounds[ir][0]) + ' ' + ' '.join(["{:5.2f}".format(b) for b in bounds[ir][1:]]) + '\n')

# Print out the bin values to be copy/pasted to input file.
print("R bins for copy/paste to 3DLIM input file:")
print("Number of bins: {}".format(len(lim_rbins)))
for r in lim_rbins:
    print("{:.3f}".format(r))
print()
print("P bins for copy/paste to 3DLIM input file:")
print("Number of bins: {}".format(len(lim_pbins)))
for p in lim_pbins:
    print("{:.3f}".format(p))
print()
nybins = 150
print("Option for Y Bins")
print("Number of bins: {}".format(nybins))
lim_ybins = np.linspace(0, bounds.max(), nybins)
for y in lim_ybins:
    print("{:.2f}".format(y))
print()
print("Some additional comments")
print("  - Change AW to something less than {:.3f}".format(lim_rbins[0]))
print("  - Change A to something greater than {:.3f}".format(lim_rbins[-1]))
print("  - Maximum bounds value is {:.3f}".format(bounds.max()))

# Read in 3D wall for below plot.
wall_path = "/Users/zamperini/Documents/d3d_work/184527/mafot/mafot_3D_wall.dat"
wall = read_wall(wall_path)

# Read in the EFIT cross section for comparison.
#efit_path = "/Users/zamperini/Documents/d3d_work/184527/SAS_prototypeL1R2b_extmesh2.ogr"
efit_path = "/Users/zamperini/Documents/d3d_work/d3d_wall_sas_fixed.dat"
efit_df = pd.read_csv(efit_path, skiprows=1, names=["R1", "Z1", "R2", "Z2"], delimiter="\t")

# A Polygon to plot showing the region we are simulating.
poly = Polygon([(min(rs), min(zs)), (min(rs), max(zs)), (max(rs), max(zs)), (max(rs), min(zs))])

# Helper plot to look at the bounds file plotted.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7,7))

ax1.plot(efit_df["R1"], efit_df["Z1"], color="k", lw=1, ls=":")
ax1.plot(wall[360-tor_origin][0], wall[360-tor_origin][1], color="k", lw=2)
ax1.plot(*poly.exterior.xy, color="r", lw=3)
ax1.set_aspect("equal")
ax1.axis("off")

vmin = 0.1; vmax = 50
levels=np.geomspace(vmin, vmax, 11)
cont = ax2.pcolormesh(rs, zs, bounds.T, shading="auto",
  norm=LogNorm(vmin, vmax, clip=True), cmap="inferno")
cbar = fig.colorbar(cont, ticks=levels, format="%3.1f")
cbar.set_label("Connection Length (m)")
ax2.set_xlabel("R (m)")
ax2.set_ylabel("Z (m)")
ax2.set_aspect("equal")
fig.tight_layout()
fig.show()
