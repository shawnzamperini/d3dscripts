import numpy as np
import postgkyl as pg
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker, cm
import time


mD = 931.49e6 / 3e8 ** 2

class Plasma3:

    def __init__(self, root, fname, fstart, fend):
        """
        root (str): Path to the directory containing the background data.
        fstart (int): Starting frame number.
        fend (end): Ending frame number.
        """

        # Load in the plasma background data. This is taken from Tess's nstx-diagnostic.py script.
        step = 1
        Nt = fend - fstart
        Nt //= step
        dt = 0.5e-6

        fp = "{:}{:}".format(root, fname)
        # fp = "{:}d3d-167196".format(root)
        print("Data location: {:}".format(fp))

        # physical constants
        mp = 1.672623e-27
        AMU = 2.014  # Deuterium ions
        mi = mp * AMU
        me = 9.10938188e-31
        eV = 1.602e-19

        c_s = np.sqrt(40 * eV / mi)

        # For equilibrium profiles
        d = pg.GData("%s_electron_M0_%d.bp" % (fp, fstart))
        dg = pg.GInterpModal(d, 1, 'ms')
        X, nElc = dg.interpolate(0)
        nElc = nElc[:, :, :, 0]

        zi = len(X[2]) // 2  # to take a slice at midplane

        Nx = len(X[0]) - 1
        Nz = len(X[2]) - 1
        Ny = len(X[1]) - 1

        # calculate grids
        x = X[0]
        dx = np.diff(x)[0]
        x = x + dx / 2
        x = x[:-1]

        y = X[1]
        dy = np.diff(y)[0]
        y = y + dy / 2
        y = y[:-1]

        z = X[2]
        dz = np.diff(z)[0]
        z = z + dz / 2
        z = z[:-1]

        # Create arrays for equilibrium profiles that are averaged in y, z and time
        nElc_tot = []
        tElc_tot = []
        nIon_tot = []
        tIon_tot = []
        nCar_tot = []
        tCar_tot = []
        phi_tot = []
        elecx_tot = []
        elecy_tot = []
        elecz_tot = []
        print("Loading background plasma...")
        for t in range(fstart, fend, step):
            if t % 100 == 0:
                print('frame = ', t)
            if (t == fend - 1):
                print("frame = ", t)

            # electron density
            d = pg.GData("%s_electron_M0_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, nElc = dg.interpolate(0)
            nElc = nElc[:, :, :, 0]
            nElc_tot.append(nElc)

            # electron temperature
            d = pg.GData("%s_electron_Temp_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, tElc = dg.interpolate(0)
            tElc = tElc[:, :, :, 0] / eV
            tElc_tot.append(tElc)

            # ion density
            d = pg.GData("%s_ion_M0_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, nIon = dg.interpolate(0)
            nIon = nIon[:, :, :, 0]
            nIon_tot.append(nIon)

            # ion temperature
            d = pg.GData("%s_ion_Temp_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, tIon = dg.interpolate(0)
            tIon = tIon[:, :, :, 0] / eV
            tIon_tot.append(tIon)

            # carbon density
            d = pg.GData("%s_carbon_M0_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, nCar = dg.interpolate(0)
            nCar = nCar[:, :, :, 0]
            nCar_tot.append(nCar)

            # carbon temperature
            d = pg.GData("%s_ion_Temp_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, tCar = dg.interpolate(0)
            tCar = tCar[:, :, :, 0] / eV
            tCar_tot.append(tCar)

            # Phi eq profile calcuation
            d = pg.GData("%s_phi_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, phi = dg.interpolate(0)
            phiZmin = phi[:, :, 0, 0]
            phi = phi[:, :, :, 0]
            phi_tot.append(phi)

            # The electric field as the gradient of the potential.
            elec = np.gradient(-phi, x, y, z)
            elecx_tot.append(elec[0])
            elecy_tot.append(elec[1])
            elecz_tot.append(elec[2])

        nElc_tot = np.array(nElc_tot)
        tElc_tot = np.array(tElc_tot)
        nIon_tot = np.array(nIon_tot)
        tIon_tot = np.array(tIon_tot)
        nCar_tot = np.array(nCar_tot)
        tCar_tot = np.array(tCar_tot)
        phi_tot = np.array(phi_tot)
        elecx_tot = np.array(elecx_tot)
        elecy_tot = np.array(elecy_tot)
        elecz_tot = np.array(elecz_tot)
        self.ne = nElc_tot
        self.te = tElc_tot
        self.ni = nIon_tot
        self.ti = tIon_tot
        self.nc = nCar_tot
        self.tc = tCar_tot
        self.vp = phi_tot
        self.er = elecx_tot
        self.ep = elecy_tot
        self.ez = elecz_tot
        self.r = x
        self.p = y
        self.z = z
        self.dt = dt
        self.cs = np.sqrt((self.te + self.ti) / mD)

        # Calculate the B field and its gradients.
        def bmag(x):
            B_axis = 2.04
            R = 2.30
            R0 = 1.722
            B0 = B_axis * (R0 / R)
            return B0 * R / x
        BB = np.zeros(self.ne.shape[1:])
        for i in range(0, len(self.r)):
            BB[i, :, :] = bmag(self.r[i])
        self.BB = BB

        print("Calculating the magnetic field gradient...")
        gradB = np.gradient(BB, x, y, z)
        self.gradBx = gradB[0]  # Units of T / m
        self.gradBy = gradB[1]  # ''
        self.gradBz = gradB[2]  # Units of T / radian


plasma = Plasma3("/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-167196-v3-c3/", "d3d-167196-v3-c3", 0, 11)

# Bounds.
rmin = plasma.r.min()
rmax = plasma.r.max()
pmin = plasma.p.min()
pmax = plasma.p.max()
midz = len(plasma.z) // 2

# Some numbers for the plots below.
ep_min = plasma.ep[:, :, :, midz].min()
ep_max = plasma.ep[:, :, :, midz].max()
ep_lim = np.max(np.abs([ep_min, ep_max])) / 2
plot_ep = np.clip(plasma.ep, -ep_lim, ep_lim)
time_label = 0.0
scat_cmap = matplotlib.cm.get_cmap("PiYG")
max_raddist = 0.03

# Plot summary. First assemble arrays of the particle position histories.
back_data = "conc_c"
conc_c = plasma.nc / plasma.ni
conc_c_lims = [conc_c.min(), conc_c.max()]
fontsize = 12
figsize = (7, 6)
fig, ax = plt.subplots(figsize=figsize)
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
if back_data == "ne":
    cont = ax.contourf(plasma.r, plasma.p, plasma.ne[0, :, :, midz].T, zorder=10)
    cbar = fig.colorbar(cont, cax=cax, extend="both")
    cbar.set_label(r"$\mathdefault{n}_{e}$ (m-3)", fontsize=fontsize)
elif back_data == "ni":
    cont = ax.contourf(plasma.r, plasma.p, plasma.ni[0, :, :, midz].T, zorder=10)
    cbar = fig.colorbar(cont, cax=cax, extend="both")
    cbar.set_label(r"$\mathdefault{n}_{i}$ (m-3)", fontsize=fontsize)
elif back_data == "nc":
    cont = ax.contourf(plasma.r, plasma.p, plasma.nc[0, :, :, midz].T, zorder=10)
    cbar = fig.colorbar(cont, cax=cax, extend="both")
    cbar.set_label(r"$\mathdefault{n}_{c}$ (m-3)", fontsize=fontsize)
elif back_data == "conc_c":
    cont = ax.contourf(plasma.r, plasma.p, conc_c[0, :, :, midz].T, zorder=10,
                       levels=np.linspace(conc_c_lims[0], conc_c_lims[1], 11))
    cbar = fig.colorbar(cont, cax=cax, extend="both")
    cbar.set_label(r"$\mathdefault{n}_{c} / \mathdefault{n}_{i}$", fontsize=fontsize)
elif back_data == "ep":
    cont = ax.contourf(plasma.r, plasma.p, plot_ep[0, :, :, midz].T, zorder=10,
                       levels=np.linspace(-ep_lim, ep_lim, 11), cmap="coolwarm")
    cbar = fig.colorbar(cont, cax=cax, extend="both")
    cbar.set_label(r"$\mathdefault{E}_{\perp}$ (V/m)", fontsize=fontsize)

ax.set_xlabel("Radial (m)", fontsize=fontsize)
ax.set_ylabel("Binormal (m)", fontsize=fontsize)
ax.set_title("{:7.2f} us".format(time_label), fontsize=fontsize)
fig.tight_layout()


def update(frame):
    # Update contour.
    ax.clear()
    if back_data == "ne":
        cont = ax.contourf(plasma.r, plasma.p, plasma.ne[frame, :, :, midz].T, zorder=10)
        cax.cla()
        cbar = fig.colorbar(cont, cax=cax, extend="both")
        cbar.set_label(r"$\mathdefault{n}_{e}$ (m-3)", fontsize=fontsize)
    elif back_data == "ni":
        cont = ax.contourf(plasma.r, plasma.p, plasma.ni[frame, :, :, midz].T, zorder=10)
        cax.cla()
        cbar = fig.colorbar(cont, cax=cax, extend="both")
        cbar.set_label(r"$\mathdefault{n}_{i}$ (m-3)", fontsize=fontsize)
    elif back_data == "nc":
        cont = ax.contourf(plasma.r, plasma.p, plasma.nc[frame, :, :, midz].T, zorder=10)
        cax.cla()
        cbar = fig.colorbar(cont, cax=cax, extend="both")
        cbar.set_label(r"$\mathdefault{n}_{c}$ (m-3)", fontsize=fontsize)
    elif back_data == "conc_c":
        cont = ax.contourf(plasma.r, plasma.p, conc_c[frame, :, :, midz].T, zorder=10,
                           levels=np.linspace(conc_c_lims[0], conc_c_lims[1], 11))
        cax.cla()
        cbar = fig.colorbar(cont, cax=cax, extend="both")
        cbar.set_label(r"$\mathdefault{n}_{c} / \mathdefault{n}_{i}$", fontsize=fontsize)
    elif back_data == "ep":
        cont = ax.contourf(plasma.r, plasma.p, plot_ep[frame, :, :, midz].T, zorder=10,
                           levels=np.linspace(-ep_lim, ep_lim, 11), cmap="coolwarm")
        cax.cla()
        cbar = fig.colorbar(cont, cax=cax, extend="both")
        cbar.set_label(r"$\mathdefault{E}_{\perp}$ (V/m)", fontsize=fontsize)


    ax.set_xlabel("Radial (m)", fontsize=fontsize)
    ax.set_ylabel("Binormal (m)", fontsize=fontsize)
    ax.set_title("frame {:}".format(frame), fontsize=fontsize)
    fig.tight_layout()
    return cont

print("Creating animation...")
start = time.time()
interval = 80  # Usual value.
# interval = 40
ani = animation.FuncAnimation(fig=fig, func=update, frames=11 - 1, interval=interval)
end = time.time()
print("Done (took {:.0f} seconds)".format(end - start))
fname = "diiid-167196-v3-c3"

start = time.time()
ani.save(filename="{:}.html".format(fname), writer="html")
end = time.time()
print("HTML saved (took {:.0f} seconds)".format(end - start))
# start = time.time()
# ani.save(filename="{:}.gif".format(fname), writer="pillow")
# end = time.time()
# print("GIF saved (took {:.0f} seconds)".format(end - start))
start = time.time()
writervideo = animation.FFMpegWriter(fps=25)
# ani.save('{:}.mp4'.format(fname), writer=writervideo)
end = time.time()
print("MP4 saved (took {:.0f} seconds)".format(end - start))
# plt.show()