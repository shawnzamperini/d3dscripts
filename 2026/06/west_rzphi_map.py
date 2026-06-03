from read_gfile import read_gfile
import matplotlib.pyplot as plt
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

def parameterize_by_theta(R_vals, Z_vals, R_axis, Z_axis, Ntheta=256):
    """
    Given contour points (R_vals, Z_vals), return R(θ), Z(θ)
    on a uniform θ grid.
    """
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

# Average over y (alpha), and since 
# phi = (theta - alpha) / q(psi) = (z - y) / q(x)
# a y-average is a toroidal average since phi ~ y at constant x, z
#nz_yavg = fp.nc["output"]["nz"][:].mean(axis=0).mean(axis=1)  # t then y average
nz_yavg = fp.nc["output"]["nz"][-1].mean(axis=1)  # last frame then y average
#nz_yavg = fp.nc["background"]["Ui_Z"][-1].mean(axis=1)  # last frame then y average

# Assemble lists of the R, Z and toroidally average nz values
Rs = []
Zs = []
nzs = []
for i in range(len(x)):
	for j in range(len(z)):
		Rs.append(R_of_psitheta((x[i], z[j])).item())
		Zs.append(Z_of_psitheta((x[i], z[j])).item())
		nzs.append(nz_yavg[i, j])
		#nzs.append(nz_yavg[xidx[i], j])
		#nzs.append(1.0)
	
# Plot of tricountourf
#fig, ax = plt.subplots()
#ax.tricontourf(Rs, Zs, nzs)
#ax.axis('equal')
#fig.show()

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
tcf = ax.tricontourf(Rs, Zs, nzs)

# Plot white over everything outside the simulation's psi range
ax.contourf(rgrid, zgrid, psirz, 
            levels=[-np.inf, grid_x.min(), grid_x.max(), np.inf],
            colors=['white', 'none', 'white'])

# Show where the cutoff is, which we call Zdiv
ax.plot([2.1, 2.6], [Zdiv, Zdiv], linestyle="--", color="k")

#ax.scatter(Rs, Zs, s=15)
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

cbar = fig.colorbar(tcf, ax=ax)
cbar.set_label("Ui_Z")

fig.show()


