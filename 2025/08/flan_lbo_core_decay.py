import flan_plots
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy.stats import linregress


def get_tot_core_imp(path, molar_mass, z_density):

	# Load LBO simulation
	fp = flan_plots.FlanPlots(path)

	# Variables needed
	t = fp.nc["time"][:] * 1e6  # s to us
	grid_x = fp.nc["x"][:]
	grid_y = fp.nc["y"][:]
	grid_z = fp.nc["z"][:]
	J = fp.nc["jacobian"][:]
	nz = fp.nc["imp_density"][:]

	# Additional settings
	x_lcfs = 0.10
	num_particles = 100000
	pellet_radius = 0.0001
	mass_pellet = 4 / 3 * np.pi * pellet_radius**3 * z_density
	num_moles = mass_pellet / molar_mass
	num_atoms = num_moles * 6.022e23

	# Loop through every cell in the core, summing up the number of total particles
	# in the core at each time frame
	tot_core_imp = np.zeros(len(t))
	for i in tqdm(range(len(t))):

		# If all zeros then nothing to be done
		if (nz[i].sum() == 0): continue
		
		for j in range(len(grid_x)-1):
			
			# Limit to just the core
			if grid_x[j] > x_lcfs: continue

			for k in range(len(grid_y)-1):
				for l in range(len(grid_z)-1):

					# Calculate volume of cell
					dx = grid_x[j+1] - grid_x[j]
					dy = grid_y[k+1] - grid_y[k]
					dz = grid_z[l+1] - grid_z[l]
					V = J[j,k,l] * dx * dy * dz

					# Calculate total number of particles (nz * V)
					tot_core_imp[i] += nz[i,j,k,l] * V

	# Need to think about how to get here, but regardless it makes sense to
	# normalize the data and multiply by the total number of atoms since they
	# all start in the core
	#tot_core_imp = tot_core_imp / tot_core_imp.max() * num_atoms
	tot_core_imp = tot_core_imp / tot_core_imp.max()

	return t, tot_core_imp

t, tot_core_he = get_tot_core_imp('/home/zamp/flandir/nt_v4/nt_v4_lbo_he.nc', 4.0026,  0.1785 * 1e3)
t, tot_core_fe = get_tot_core_imp('/home/zamp/flandir/nt_v4/nt_v4_lbo_fe.nc', 55.845,  7874 * 1e3)
t, tot_core_w  = get_tot_core_imp('/home/zamp/flandir/nt_v4/nt_v4_lbo_w.nc',  183.840, 19250 * 1e3)

# Wrap up into dictionary
tot_core_imp = {"He":tot_core_he, "Fe":tot_core_fe, "W":tot_core_w}

# Exponential fit to get 1/e decay length for indicated region
def get_decay_length(x, xmin, xmax, y):
	
	# Mask data
	mask = np.logical_and(x >= xmin, x <= xmax)
	x = x[mask]
	y = y[mask]

	# Return error if zero in y data
	if (0 in y): print("Error! Can't fit if y contains zeros")

	# Natural log of data
	ln_y = np.log(y)
	
	# Linear regression ln(y) = bx + ln(a)
	slope, intercept, r_value, p_value, std_err = linregress(x, ln_y)

	# Go back to exponential
	a_fit = np.exp(intercept)	
	b_fit = slope

	# Return decay length and fit values
	return 1/b_fit, x, a_fit * np.exp(b_fit * x)

# Get fit data in selected range
fit_tmin = 470
fit_tmax = 520
decay_time_he, t_fit, tot_fit_he = get_decay_length(t, fit_tmin, fit_tmax, tot_core_he)
decay_time_fe, t_fit, tot_fit_fe = get_decay_length(t, fit_tmin, fit_tmax, tot_core_fe)
decay_time_w,  t_fit, tot_fit_w  = get_decay_length(t, fit_tmin, fit_tmax, tot_core_w)

# Wrap up into dictionary
decay_times = {"He":decay_time_he, "Fe":decay_time_fe, "W":decay_time_w}

decay_labels = {}
for z in tot_core_imp.keys():
	decay_labels[z] = (fr"$\mathdefault{{\tau_{{{z}}}}}$ = {int(-decay_times[z])} $\mathdefault{{\mu s}}$")

# List of colors to cycle through
color_count = 0
colors = ["tab:red", "tab:purple", "tab:cyan", "tab:green", "tab:pink"]

# A nice plot
fontsize = 18
lw = 3
fig, ax = plt.subplots(figsize=(7, 6))

for z in tot_core_imp.keys():
	ax.plot(t-t.min(), tot_core_imp[z], lw=lw+1, color="k")
	ax.plot(t-t.min(), tot_core_imp[z], lw=lw, color=colors[color_count], label=decay_labels[z])
	color_count += 1
#ax.plot(t_fit-t.min(), tot_fit, color="tab:purple", linestyle="--", lw=lw)
ax.set_xlabel("Time (us)", fontsize=fontsize)
ax.set_ylabel("Atoms in core (norm.)", fontsize=fontsize)
ax.tick_params(labelsize=fontsize-1, which="both")
#ax.text(0.6, 0.6, decay_label_fe, transform=ax.transAxes, fontsize=fontsize)
ax.legend(fontsize=fontsize)
fig.tight_layout()
fig.show()
