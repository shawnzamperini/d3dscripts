import flan_plots
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


# Some constants.
amu = 1.660e-27  # kg
mi = amu * 2.014  # kg
mW = amu * 183.38 
me = 9.10938188e-31  # kg
ev = 1.602e-19  # C
lnalpha = 15
eps0 = 8.85e-12

ZW = 1  # W charge state
Z = 1  # Deuterium
Ti = 10
ni = 1e19
B = 1.5  # T

# Load FlanPlots object.
fp_path = "/Users/zamperini/flandir/reg_testcase1/saved_results_5/reg_testcase1.nc"
fp = flan_plots.FlanPlots(fp_path)

# Then we want the density for the plot, and the ion temparature
# for the gyro-orbits.
f = 50
ne_data = fp.plot_profiles(["ne"], plot_z=0.3125, normtype=["log"], 
    vmin=[1e18], vmax=[2e19], xlabel=r"$\mathdefault{R-R_{sep}}$ (m)", 
    ylabel="Binormal (m)", x_offset=-2.259, f=f)
ti_data = fp.plot_profiles(["ti"], plot_z=0.3125, normtype=["log"], 
    xlabel=r"$\mathdefault{R-R_{sep}}$ (m)", 
    ylabel="Binormal (m)", x_offset=-2.259, f=f)
ex_data = fp.plot_profiles(["ex"], plot_z=0.3125, normtype=["log"], 
    xlabel=r"$\mathdefault{R-R_{sep}}$ (m)", 
    ylabel="Binormal (m)", x_offset=-2.259, f=f)

# Now a plot of the electron density, like usual. 
data_f = ne_data["data"][0][f].T
data_f_ex = ex_data["data"][0][f].T
data_f_ti = ti_data["data"][0][f].T
x = ne_data["x"]
y = ne_data["y"]

# Calculate a few gyro-orbit circles at arbitrary locations.
def gyro_orbit_patch(x_gyro, y_gyro):
    
    # Find nearest indices first. 
    xidx = np.argmin(np.abs(x_gyro - x))
    yidx = np.argmin(np.abs(y_gyro - y))

    # Get Ti value and calculate gyroradius.
    ti_gyro = data_f_ti[xidx, yidx]
    vtW = np.sqrt(2 * ti_gyro * ev / mW)
    rhoW = vtW / (ZW * ev * B / mW)

    # Now create a circle patch object of this gyroorbit. 
    patch = mpl.patches.Circle((x_gyro, y_gyro), radius=rhoW, 
        edgecolor="tab:green", fill=False, lw=3)
    return patch


# Create normalization based on input options.
#norm = mpl.colors.LogNorm(vmin=1e18, vmax=2e19)
#cmap = "inferno"
norm = mpl.colors.Normalize(vmin=-5000, vmax=5000)
cmap = "coolwarm"


# Make a 2D plot using pcolormesh.
fontsize = 18
fig, ax = plt.subplots(figsize=(7, 6))
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
cont = ax.pcolormesh(x, y, data_f_ex, shading="nearest", 
    norm=norm, cmap=cmap)
ax.set_xlabel(r"$\mathdefault{R-R_{sep}}$ (m)", fontsize=fontsize)
ax.set_ylabel("Binormal (m)", fontsize=fontsize)
#ax.set_title("t = {:10.2f} us".format(dt * f * 1e6))
ax.tick_params(axis='both', which='major', labelsize=fontsize-2)

# Colorbar.
#div = make_axes_locatable(axs[i])
#cax = div.append_axes('right', '5%', '5%')
cbar = fig.colorbar(cont, cax=cax)
#cbar.set_label(r"$\mathdefault{n_e\ (m^{-3}}$)", fontsize=fontsize)
cbar.set_label(r"$\mathdefault{E_r\ (V/m)}$", fontsize=fontsize)

# Add gyro-orbit patch.
ax.add_patch(gyro_orbit_patch(0.07, 0))
ax.set_aspect("equal")

fig.tight_layout()
fig.show()
