import flan_plots
import matplotlib.pyplot as plt
import postgkyl as pg
import numpy as np


# Load FlanPlots objects
ncpath_w = "/mnt/c/Users/Shawn Zamperini/Documents/perlmutter_staging/west_antenna/antenna_coll_W.nc"
ncpath_b = "/mnt/c/Users/Shawn Zamperini/Documents/perlmutter_staging/west_antenna/antenna_coll_B.nc"
fp_w = flan_plots.FlanPlots(ncpath_w)
fp_b = flan_plots.FlanPlots(ncpath_b)

# Some arrays 
x = fp_w.nc["geometry"]["x"][:]
y = fp_w.nc["geometry"]["y"][:]
z = fp_w.nc["geometry"]["z"][:]

# The following code has a lot copied from FlanPlots. The purpose is to get
# average R, Z values for the grid so that we can parametriz the distance on 
# the poloidal plane for plots

# Toroidal average
data_yavg_w = fp_w.load_data_frame("nz", -1).mean(axis=1)
data_yavg_b = fp_b.load_data_frame("nz", -1).mean(axis=1)
dr_yavg_w = fp_w.load_data_frame("Dr", -1).mean(axis=1)
dr_yavg_b = fp_b.load_data_frame("Dr", -1).mean(axis=1)

# Load nodes file
nodes_path = "/mnt/c/Users/Shawn Zamperini/Documents/perlmutter_staging/west_antenna/gk_west_lsn_sol_3x2v_p1-nodes-hires.gkyl"
#nodes_path = "/mnt/c/Users/Shawn Zamperini/Documents/perlmutter_staging/west_antenna/gk_west_lsn_sol_3x2v_p1-nodes.gkyl"
nodes_gdata = pg.GData(nodes_path)
nodes = nodes_gdata.get_values()

# R, Z, phi coordinates stored in last dimension
R_nodes = nodes[:, :, :, 0].mean(axis=1)
Z_nodes = nodes[:, :, :, 1].mean(axis=1)
#phi_nodes = nodes[:, :, :, 2].mean(axis=1)

# Get the grid this nodes file used (which does not have to be the
# grid the simulation used)
x_nodes = nodes_gdata.get_grid()[0]
z_nodes = nodes_gdata.get_grid()[2]
x_centers = 0.5 * (x_nodes[:-1] + x_nodes[1:])
z_centers = 0.5 * (z_nodes[:-1] + z_nodes[1:])

Nx = len(x_nodes) - 1
Nz = len(z_nodes) - 1

R_node_grid = nodes[:, 0, :, 0]
Z_node_grid = nodes[:, 0, :, 1]

# Needed to like extend to the outermost cells in each dimesnion or
# something
def extend_node_grid(G):
	last_row = 2 * G[-1, :] - G[-2, :]
	G = np.vstack([G, last_row[np.newaxis, :]])
	last_col = 2 * G[:, -1] - G[:, -2]
	G = np.hstack([G, last_col[:, np.newaxis]])
	return G

R_ext = extend_node_grid(R_node_grid)
Z_ext = extend_node_grid(Z_node_grid)

# Create vertices of each cell that are used to make the polygons later
# Essentially making an irregular grid one cell at a time.
avg_R = np.zeros((len(x_centers), len(z_centers)))
avg_Z = np.zeros((len(x_centers), len(z_centers)))
nw_rz = np.zeros((len(x_centers), len(z_centers)))
nb_rz = np.zeros((len(x_centers), len(z_centers)))
Dr_w_rz = np.zeros((len(x_centers), len(z_centers)))
Dr_b_rz = np.zeros((len(x_centers), len(z_centers)))
s_pol = np.zeros((len(x_centers), len(z_centers)))
for i in range(Nx):
	for j in range(Nz - 1):

		# Average R, Z coordinates
		avg_R[i, j] = (R_ext[i,j] + R_ext[i,j+1] + R_ext[i+1,j+1] + R_ext[i+1, j]) / 4.0
		avg_Z[i, j] = (Z_ext[i,j] + Z_ext[i,j+1] + Z_ext[i+1,j+1] + Z_ext[i+1, j]) / 4.0

		# Grab the nearest Flan simulation value to this cell.
		i_sim = np.argmin(np.abs(x - x_centers[i]))
		j_sim = np.argmin(np.abs(z - z_centers[j]))
		nw_rz[i, j] = data_yavg_w[i_sim, j_sim]
		nb_rz[i, j] = data_yavg_b[i_sim, j_sim]
		Dr_w_rz[i, j] = dr_yavg_w[i_sim, j_sim]
		Dr_b_rz[i, j] = dr_yavg_b[i_sim, j_sim]

		# Parametrize the data from the first z coordinate to the last
		if j != 0:
			dist = np.sqrt(np.square(avg_R[i,j-1] - avg_R[i,j]) 
				+ np.square(avg_Z[i,j-1] - avg_Z[i,j]))
			s_pol[i, j] = s_pol[i, j-1] + dist


# Separatrix density vs poloidal s coordinate. Ditch last point since it's
# a (0,0), extra index snuck in there.
xidx = -1  # -1 = separatrix, 0 = edge
avg_R_sep = avg_R[xidx, :-1]
avg_Z_sep = avg_Z[xidx, :-1]
nw_pol_sep = nw_rz[xidx, :-1]
nb_pol_sep = nb_rz[xidx, :-1]
Dr_w_pol_sep = Dr_w_rz[xidx, :-1]
Dr_b_pol_sep = Dr_b_rz[xidx, :-1]
s_pol_sep = s_pol[xidx, :-1]

# Locate the OMP and IMP
Rmaxis = 2.547431773
Zmaxis = 0.03770508942
out_mask = avg_R_sep > Rmaxis
in_mask = avg_R_sep < Rmaxis
omp_idx = np.argmin(np.abs(avg_Z_sep[out_mask] - Zmaxis))
imp_idx = np.argmin(np.abs(avg_Z_sep[in_mask] - Zmaxis))
s_pol_sep_omp = s_pol_sep[out_mask][omp_idx]
s_pol_sep_imp = s_pol_sep[in_mask][imp_idx]

# Get a psin value
psi_sep = 0.4167022815
psi_axis = 0.5777173699
psi_norm = (x_centers[xidx] - psi_axis) / (psi_sep - psi_axis)

# Impurity density
fontsize = 16
fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.plot(s_pol_sep, nw_pol_sep, label="W")
ax1.plot(s_pol_sep, nb_pol_sep, label="B")
ax1.axvline(s_pol_sep_omp, color="k", linestyle="--")
ax1.axvline(s_pol_sep_imp, color="r", linestyle="--")
ax1.set_xlabel("Poloidal distance from outer target (m)", fontsize=fontsize)
ax1.set_ylabel("nZ (arb.)", fontsize=fontsize)
ax1.set_title("psi_n = {:.3f}".format(psi_norm), fontsize=fontsize)
ax1.legend(fontsize=fontsize)
fig.tight_layout()
fig.show()

# Radial diffusion coefficient
fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.plot(s_pol_sep, Dr_w_pol_sep, label="W")
ax1.plot(s_pol_sep, Dr_b_pol_sep, label="B")
ax1.axvline(s_pol_sep_omp, color="k", linestyle="--")
ax1.axvline(s_pol_sep_imp, color="r", linestyle="--")
ax1.set_xlabel("Poloidal distance from outer target (m)", fontsize=fontsize)
ax1.set_ylabel("Dr (m2/s)", fontsize=fontsize)
ax1.set_title("psi_n = {:.3f}".format(psi_norm), fontsize=fontsize)
ax1.legend(fontsize=fontsize)

fig.tight_layout()
fig.show()

#fig, ax1 = plt.subplots(figsize=(5, 4))
#ax1.plot(s_pol_sep, avg_R_sep, label="R")
#ax1.plot(s_pol_sep, avg_Z_sep, label="Z")
#ax1.set_xlabel("Poloidal distance from outer target (m)", fontsize=fontsize)
#ax1.set_title("psi = {:.3f}".format(psi_norm), fontsize=fontsize)
#ax1.legend(fontsize=fontsize)
#fig.tight_layout()
#fig.show()
