from read_gfile import read_gfile
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import flan_plots
from skimage import measure
from scipy.interpolate import interp1d, RegularGridInterpolator, griddata


# First load gfile for WEST #62104 and pull out the R, Z grid and the 2D psi
# array defined on it
path = "west_LSN.geqdsk"
gfile = read_gfile(path)
rgrid = gfile["rgrid"]
zgrid = gfile["zgrid"]
psigrid = gfile["psi_grid"]
psirz = gfile["psirz"].T  # Indexing Fortran --> C (which is assumed in numpy)
R_axis = gfile["rmaxis"]
Z_axis = gfile["zmaxis"]
psi_axis = gfile["simag"]
psi_lcfs = gfile["sibry"]
qpsi = gfile["qpsi"]
rlim = gfile["rlim"]
zlim = gfile["zlim"]

# To simplify the geometry, Tess set the floor at a constant Z = -0.6 m. So in
# our plot we should remove all data below this. It ends up causing issues
# to have data below this because it lands data points on a countour outside
# the vessel and simulation volume. 
Zdiv = -0.6
Zmask = zgrid > Zdiv
zgrid = zgrid[Zmask]
psirz = psirz[Zmask, :]

# Then load Flan simulation
#path = "/mnt/c/Users/Shawn Zamperini/Documents/gkyldir_staging/ot_point_source.nc"
path = "/pscratch/sd/z/zamp/flandir/west_nearsol_point/ot_point_source.nc"
fp = flan_plots.FlanPlots(path)

# Pull out the grid coordinates
x = fp.nc["geometry"]["x"][:]  # psi
y = fp.nc["geometry"]["y"][:]  # alpha
z = fp.nc["geometry"]["z"][:]  # chi or theta
grid_x = fp.nc["geometry"]["grid_x"][:]  # psi

# Limit to a spcific x range
#xidx = [-1]  # 0 = SOL edge, -1 = separatrix
#x = x[xidx]

# ------------------------------------------------------------
# 2. Extract a single flux surface ψ = ψ0
# ------------------------------------------------------------

def extract_flux_surface(psirz, R, Z, psi0):
    """
    Returns R_contour, Z_contour for the flux surface ψ = ψ0.
    """
    # skimage.find_contours expects array indexed as [Z, R]
	# Each contour is in index space, it doesn't know about the R, Z
	# coordinates. Fractional indices can be returned! So if you see 
	# Ridx = 34.7, it means you should interpolate the values between idx
	# 34 and 35 to see what value 34.7 gives you.
    contours = measure.find_contours(psirz, psi0)

    if len(contours) == 0:
        raise RuntimeError(f"No contour found for psi={psi0}")

    # Choose the longest contour (usually the closed flux surface)
    contour = max(contours, key=len)

    # Convert from index space to physical coordinates
    Z_idx = contour[:, 0]
    R_idx = contour[:, 1]
    #R_idx = contour[:, 0]
    #Z_idx = contour[:, 1]

    R_vals = np.interp(R_idx, np.arange(len(R)), R)
    Z_vals = np.interp(Z_idx, np.arange(len(Z)), Z)

    return R_vals, Z_vals

# ------------------------------------------------------------
# 3. Parameterize a flux surface by θ
# ------------------------------------------------------------

"""
def parameterize_by_theta(R_vals, Z_vals, R_axis, Z_axis, Ntheta=256):
    #Given contour points (R_vals, Z_vals), return R(θ), Z(θ)
    #on a uniform θ grid.
    theta = np.arctan2(Z_vals - Z_axis, R_vals - R_axis)

    # Sort by θ
    idx = np.argsort(theta)
    theta_sorted = theta[idx]
    R_sorted = R_vals[idx]
    Z_sorted = Z_vals[idx]

    # Uniform θ grid
    theta_grid = np.linspace(-np.pi, np.pi, Ntheta)

    # Interpolate
    R_theta = interp1d(theta_sorted, R_sorted, fill_value="extrapolate")(theta_grid)
    Z_theta = interp1d(theta_sorted, Z_sorted, fill_value="extrapolate")(theta_grid)

    return theta_grid, R_theta, Z_theta
"""

def parameterize_by_theta(R_vals, Z_vals, R_axis, Z_axis, Ntheta=256):
    theta = np.arctan2(Z_vals - Z_axis, R_vals - R_axis)

    idx = np.argsort(theta)
    theta_sorted = theta[idx]
    R_sorted = R_vals[idx]
    Z_sorted = Z_vals[idx]

    # Close the periodic gap by appending the first point shifted by 2*pi
    #theta_sorted = np.append(theta_sorted, theta_sorted[0] + 2*np.pi)
    #R_sorted = np.append(R_sorted, R_sorted[0])
    #Z_sorted = np.append(Z_sorted, Z_sorted[0])

    # Wrap both ends to guarantee coverage of [-pi, pi]
    theta_sorted = np.concatenate([[theta_sorted[0] - 2*np.pi], theta_sorted, [theta_sorted[-1] + 2*np.pi]])
    R_sorted = np.concatenate([[R_sorted[-1]], R_sorted, [R_sorted[0]]])
    Z_sorted = np.concatenate([[Z_sorted[-1]], Z_sorted, [Z_sorted[0]]])

    theta_grid = np.linspace(-np.pi, np.pi, Ntheta)

    R_theta = interp1d(theta_sorted, R_sorted)(theta_grid)
    Z_theta = interp1d(theta_sorted, Z_sorted)(theta_grid)

    return theta_grid, R_theta, Z_theta

# ------------------------------------------------------------
# 4. Build R(ψ,θ) and Z(ψ,θ)
# ------------------------------------------------------------

def build_flux_coordinate_mapping(filename, Npsi=64, Ntheta=256):
    g, R, Z, psirz, R_axis, Z_axis, psi_axis, psi_lcfs = load_eq(filename)

    # Choose ψ grid
    psi_grid = np.linspace(psi_axis, psi_lcfs, Npsi)

    # Allocate arrays
    R_psitheta = np.zeros((Npsi, Ntheta))
    Z_psitheta = np.zeros((Npsi, Ntheta))

    # Loop over ψ surfaces
    for i, psi0 in enumerate(psi_grid):
        R_vals, Z_vals = extract_flux_surface(psirz, R, Z, psi0)
        theta_grid, R_theta, Z_theta = parameterize_by_theta(
            R_vals, Z_vals, R_axis, Z_axis, Ntheta
        )
        R_psitheta[i, :] = R_theta
        Z_psitheta[i, :] = Z_theta

    # Build interpolators
    R_of_psitheta = RegularGridInterpolator(
        (psi_grid, theta_grid), R_psitheta, bounds_error=False, fill_value=None
    )
    Z_of_psitheta = RegularGridInterpolator(
        (psi_grid, theta_grid), Z_psitheta, bounds_error=False, fill_value=None
    )

    return psi_grid, theta_grid, R_psitheta, Z_psitheta, R_of_psitheta, Z_of_psitheta


# Allocate arrays
R_psitheta = np.zeros((len(x), len(z)))
Z_psitheta = np.zeros((len(x), len(z)))

# Loop through one psi value at a time
for i, psi in enumerate(x):

	# Pull out R, Z contours of each flux tube, defined by psi value
	R_vals, Z_vals = extract_flux_surface(psirz, rgrid, zgrid, psi)

	# Parameterize the R and Z values of the flux surface by theta
	theta_grid, R_theta, Z_theta = parameterize_by_theta(R_vals, Z_vals, R_axis, Z_axis, len(z))

	R_psitheta[i, :] = R_theta
	Z_psitheta[i, :] = Z_theta

# Build interpolators
R_of_psitheta = RegularGridInterpolator(
	(x, z), R_psitheta, bounds_error=False, fill_value=None)
Z_of_psitheta = RegularGridInterpolator(
	(x, z), Z_psitheta, bounds_error=False, fill_value=None)

# ---------------------------------
import postgkyl as pg

# Important to load the high resolution nodes!
data = pg.GData("/global/cfs/cdirs/m3739/gkeyll/gkyl_for_flan/west-sol-62104/Nz128/gk_west_lsn_sol_3x2v_p1-nodes.gkyl")
data
print(data.get_grid())
print(data.get_values().shape)

nodes = data.get_values()  # shape (17, 9, 25, 3)

# nodes[i, :, k, :] gives the physical coords at all y nodes
# for psi index i and z index k.
# Average over y (axis=1) to get the R, Z at each (psi, z) point
R_nodes = nodes[:, :, :, 0].mean(axis=1)  # shape (17, 25)
Z_nodes = nodes[:, :, :, 1].mean(axis=1)  # shape (17, 25)

# x and z node coordinates
x_nodes = data.get_grid()[0]  # length 17 - actually 18
z_nodes = data.get_grid()[2]  # length 25 - actually 26
x_centers = 0.5 * (x_nodes[:-1] + x_nodes[1:])  # length 17
z_centers = 0.5 * (z_nodes[:-1] + z_nodes[1:])  # length 25

from scipy.interpolate import RegularGridInterpolator
R_of_xz = RegularGridInterpolator((x_centers, z_centers), R_nodes, 
	bounds_error=False, fill_value=None)
Z_of_xz = RegularGridInterpolator((x_centers, z_centers), Z_nodes, 
	bounds_error=False, fill_value=None)



# Build interpolators
#R_of_xz = RegularGridInterpolator((x_nodes, z_nodes), R_nodes)
#Z_of_xz = RegularGridInterpolator((x_nodes, z_nodes), Z_nodes)

# Average over y (alpha), and since 
# phi = (theta - alpha) / q(psi) = (z - y) / q(x)
# a y-average is a toroidal average since phi ~ y at constant x, z
#nz_yavg = fp.nc["output"]["nz"][:].mean(axis=0).mean(axis=1)  # t then y average
#nz_yavg = fp.nc["output"]["nz"][-1].mean(axis=1)  # last frame then y average
nz_yavg = fp.nc["output"]["nz"][0].mean(axis=1)  # first frame then y average
#nz_yavg = fp.nc["background"]["Ui_Z"][-1].mean(axis=1)  # last frame then y average

# Assemble lists of the R, Z and toroidally average nz values
Rs = []
Zs = []
nzs = []
xs = []
zs = []
for i in range(len(x)):
	for j in range(len(z)):
		#Rs.append(R_of_psitheta((x[i], z[j])).item())
		#Zs.append(Z_of_psitheta((x[i], z[j])).item())
		Rs.append(R_of_xz((x[i], z[j])).item())
		Zs.append(Z_of_xz((x[i], z[j])).item())
		nzs.append(nz_yavg[i, j])
		#nzs.append(nz_yavg[xidx[i], j])
		#nzs.append(1.0)
		xs.append(x[i])
		zs.append(z[j])

# Again, but for the node coordinates
Rs_nodes = []
Zs_nodes = []
for i in range(len(x_nodes)):
	for j in range(len(z_nodes)):
		Rs_nodes.append(R_of_xz((x_nodes[i], z_nodes[j])).item())
		Zs_nodes.append(Z_of_xz((x_nodes[i], z_nodes[j])).item())
	
# set zeros to nan's so the logscale doesn't get tripped up by the zeros
nzs = np.where(np.array(nzs) <= 0, np.nan, nzs)

mask = np.isfinite(nzs)
Rs2 = np.array(Rs)[mask]
Zs2 = np.array(Zs)[mask]
nzs2 = np.array(nzs)[mask]
xs2 = np.array(xs)[mask]
zs2 = np.array(zs)[mask]

# Plot of tricountourf
#fig, ax = plt.subplots()
#ax.tricontourf(Rs, Zs, nzs)
#ax.axis('equal')
#fig.show()

# Starting location
start_R = R_of_xz((0.4095, -3.1415))
start_Z = Z_of_xz((0.4095, -3.1415))


# Create a mesh

# 1. Make a regular grid
Rgrid = np.linspace(min(Rs), max(Rs), 500)
Zgrid = np.linspace(min(Zs), max(Zs), 500)
RR, ZZ = np.meshgrid(Rgrid, Zgrid)

# 2. Interpolate nz onto this grid
NZgrid = griddata((Rs, Zs), nzs, (RR, ZZ), method='linear')

# 3. Plot
fig, ax = plt.subplots()
#ax.pcolormesh(RR, ZZ, NZgrid, shading='auto')
#tcf = ax.tricontourf(Rs, Zs, nzs, levels=np.linspace(-10000, 10000), cmap="coolwarm")
#tcf = ax.tricontourf(Rs, Zs, nzs)
#tcf = ax.tricontourf(
#    Rs2, Zs2, nzs2,
#    norm=LogNorm(vmin=np.nanmin(nzs), vmax=np.nanmax(nzs)),
#    levels=np.logspace(np.log10(np.nanmin(nzs)), np.log10(np.nanmax(nzs)), 50),
#    cmap="inferno"
#)

#ax.scatter(start_R, start_Z, s=75, marker="*")


# Plot white over everything outside the simulation's psi range
#ax.contourf(rgrid, zgrid, psirz, 
#            levels=[-np.inf, grid_x.min(), grid_x.max(), np.inf],
#            colors=['white', 'none', 'white'])

# Show where the cutoff is, which we call Zdiv
ax.plot([2.1, 2.6], [Zdiv, Zdiv], linestyle="--", color="k")

# Labels of (x,z). Only visible if you zoom in
#ax.scatter(Rs, Zs, s=15)
#for Ri, Zi, xi, zi in zip(Rs2, Zs2, xs2, zs2):
#    ax.text(Ri, Zi, f"({xi:.4f}, {zi:.2f})",
#            fontsize=8, color="k",
#            ha="center", va="center", zorder=1e10)

ax.contour(rgrid, zgrid, psirz, levels=[psi_lcfs], colors="r")
ax.plot(rlim, zlim, color="k")
ax.axis('equal')

# Plot of the poloidal angle to help when referring where to start particles
nang = 8
line_len = 1.0
for i in range(nang):
	angle = 2.0 * np.pi * i / nang  

	# End point of line
	R1 = R_axis + line_len * np.cos(angle)
	Z1 = Z_axis + line_len * np.sin(angle)

	ax.plot([R_axis, R1], [Z_axis, Z1], color="k", linestyle="--")

#cbar = fig.colorbar(tcf, ax=ax)
#cbar.set_label("nZ")

fig.show()

# After building R_psitheta and Z_psitheta, check what R,Z
# a point at z = pi/2 maps to for your starting psi value
i_start = 0  # whatever psi index your particles start at
j_start = np.argmin(np.abs(z - np.pi/2))

print(f"z[j_start] = {z[j_start]}")
print(f"R at (psi_start, pi/2) = {R_psitheta[i_start, j_start]}")
print(f"Z at (psi_start, pi/2) = {Z_psitheta[i_start, j_start]}")

#fig, ax = plt.subplots(
#ax.plot(R_psitheta[i_start, :], Z_psitheta[i_start, :], 'b-')
#ax.contour(rgrid, zgrid, psirz, levels=[x[i_start]], colors='r')
#ax.plot(rlim, zlim, color='k')
#ax.axis('equal')
#fig.show()


# --- START polygon plot ---

from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
import numpy as np

# Reshape node R, Z into 2D arrays: shape (Nx+1, Nz+1)
Nx = len(x_nodes) - 1  # number of cells in x
Nz = len(z_nodes) - 1  # number of cells in z

# R_nodes and Z_nodes are already shape (Nx, Nz) = cell centers from y-averaging
# But we need node positions at the node grid, not cell centers.
# Use the raw nodes array directly (before averaging over y):
# nodes shape is (Nx, Ny, Nz, 3), take y=0 (all y are identical as we confirmed)
R_node_grid = nodes[:, 0, :, 0]  # shape (Nx, Nz) - but these are NODE positions
Z_node_grid = nodes[:, 0, :, 1]  # shape (Nx, Nz)

# nodes array has shape matching x_nodes[:-1] x z_nodes[:-1]?
# Check:
print(f"nodes shape: {nodes.shape}")         # e.g. (17, 9, 25, 3)
print(f"x_nodes length: {len(x_nodes)}")     # e.g. 18
print(f"z_nodes length: {len(z_nodes)}")     # e.g. 26
# So nodes gives R,Z at the LOW corner of each cell, not at all 4 corners.
# We need to construct the 4 corners from adjacent cells.

# Build quads: for each cell (i,j), the 4 corners in (R,Z) are:
# (i,j), (i+1,j), (i+1,j+1), (i,j+1)
# We need R,Z at all node indices including the last row/col.
# Extend R_node_grid to (Nx+1, Nz+1) by extrapolating the boundary:

def extend_node_grid(G):
    """Extend a (Nx, Nz) node grid to (Nx+1, Nz+1) by linear extrapolation."""
    # Add one row at the end (x direction)
    last_row = 2 * G[-1, :] - G[-2, :]
    G = np.vstack([G, last_row[np.newaxis, :]])
    # Add one column at the end (z direction)
    last_col = 2 * G[:, -1] - G[:, -2]
    G = np.hstack([G, last_col[:, np.newaxis]])
    return G

R_ext = extend_node_grid(R_node_grid)  # shape (Nx+1, Nz+1)
Z_ext = extend_node_grid(Z_node_grid)  # shape (Nx+1, Nz+1)

# Now build the polygon list
verts = []
colors = []

for i in range(Nx-1):
    for j in range(Nz-1):
        # 4 corners, ordered as a closed polygon
        quad = np.array([
            [R_ext[i,   j  ], Z_ext[i,   j  ]],
            [R_ext[i+1, j  ], Z_ext[i+1, j  ]],
            [R_ext[i+1, j+1], Z_ext[i+1, j+1]],
            [R_ext[i,   j+1], Z_ext[i,   j+1]],
        ])
        verts.append(quad)
        #colors.append(nz_yavg[i, j])
        colors.append(1.0)

colors = np.array(colors)

# Mask out non-positive values for log scale
valid = colors > 0
verts_valid = [v for v, ok in zip(verts, valid) if ok]
colors_valid = colors[valid]

# Plot
fig, ax = plt.subplots()

norm = LogNorm(vmin=np.nanmin(colors_valid), vmax=np.nanmax(colors_valid))
coll = PolyCollection(verts_valid, array=colors_valid, cmap="inferno",
                      norm=norm, edgecolors='k')
ax.add_collection(coll)

ax.contour(rgrid, zgrid, psirz, levels=[psi_lcfs], colors="r")
ax.plot(rlim, zlim, color="k")
#ax.contourf(rgrid, zgrid, psirz,
#            levels=[-np.inf, grid_x.min(), grid_x.max(), np.inf],
#            colors=['white', 'none', 'white'])
ax.plot([2.1, 2.6], [Zdiv, Zdiv], linestyle="--", color="k")
ax.axis('equal')
ax.autoscale_view()

cbar = fig.colorbar(coll, ax=ax)
cbar.set_label("nZ")
fig.show()

# --- END polygon plot ---
