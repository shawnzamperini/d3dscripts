import oedge_plots
import numpy as np
import matplotlib.pyplot as plt


# Some constants.
ring = 40
mz_amu = 183.84
md_amu = 2.01
mz = mz_amu * 1.66e-27  # Mass of W in kg
md = md_amu * 1.66e-27   # Mass of D in kg
col_log = 15

# Paths to cases.
ncpath00 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-025a.nc'  # No flows.
ncpath01 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-026b.nc'  # M = -0.1
ncpath02 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-026c.nc'  # M = -0.2
ncpath03 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-026d.nc'  # M = -0.3
ncpath04 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-026e.nc'  # M = -0.4
op00 = oedge_plots.OedgePlots(ncpath00)
op01 = oedge_plots.OedgePlots(ncpath01)
op02 = oedge_plots.OedgePlots(ncpath02)
op03 = oedge_plots.OedgePlots(ncpath03)
op04 = oedge_plots.OedgePlots(ncpath04)

# Along ring quantities that don't depend on the Mach number option.
s, ti  = op00.along_ring(ring, 'KTIBS', plot_it=False)
s, ne  = op00.along_ring(ring, 'KNBS',  plot_it=False)
charges = np.arange(1, 30+1)
#charges = [10]

# The parallel background velocities (vi) for each case.
s, vi00  = op00.along_ring(ring, 'Velocity', plot_it=False)
s, vi01  = op01.along_ring(ring, 'Velocity', plot_it=False)
s, vi02  = op02.along_ring(ring, 'Velocity', plot_it=False)
s, vi03  = op03.along_ring(ring, 'Velocity', plot_it=False)
s, vi04  = op04.along_ring(ring, 'Velocity', plot_it=False)

# Array to hold the final velavg for each case.
velavg_finals = np.zeros((5, len(s)))

# Loop for each op case.
count = 0
for op in [op00, op01, op02, op03, op04]:

    # Arrays to hold results for each charge state.
    vti_all    = np.zeros((len(charges), len(s)))
    velavg_all = np.zeros((len(charges), len(s)))
    gamma_all  = np.zeros((len(charges), len(s)))
    nz_all     = np.zeros((len(charges), len(s)))

    for i in range(0, len(charges)):
        charge = charges[i]

        # Get forces.
        #s, ff = op.along_ring(ring, 'ff', plot_it=False, charge=charge)
        s, fig = op.along_ring(ring, 'fig',   plot_it=False, charge=charge)

        # Slowing down time.
        tau_s = 1.47E13 * mz_amu * ti * np.sqrt(ti / md_amu) / \
                ((1 + md_amu / mz_amu) * ne * np.power(charge, 2) * col_log)

        # vTi is FiG * tau_s / mz.
        vti = fig * tau_s / mz

        # The average velocity for this W charge state on this ring.
        s, velavg = op.along_ring(ring, 'VELAVG', plot_it=False, charge=charge)
        s, nz     = op.along_ring(ring, 'DDLIMS', plot_it=False, charge=charge)
        velavg_all[i] = velavg
        nz_all[i]     = nz
        gamma_all[i]  = velavg * nz
        vti_all[i]    = vti

    # Get average, ignoring zeros. Don't need to do for vti.
    velavg_all[velavg_all==0] = np.nan
    nz_all[nz_all==0] = np.nan
    gamma_all[gamma_all==0] = np.nan
    velavg_final = np.nanmean(velavg_all, axis=0)
    nz_final = np.nanmean(nz_all, axis=0)
    gamma_final = np.nanmean(gamma_all, axis=0)

    # Store for this case.
    velavg_finals[count] = velavg_final
    count += 1

# This won't change for any of the cases so do it after everyhting.
vti_final = vti_all.mean(axis=0)

# Index just for shorter variable names. Fix when more cases are added.
vz00 = velavg_finals[0]
vz01 = velavg_finals[1]
vz02 = velavg_finals[2]
vz03 = velavg_finals[3]
vz04 = velavg_finals[4]

# Get R-Rsep OMP.
rminrsep_omp = op00.nc.variables["MIDIST"][:].data[1][ring] * 1000
text_str = "R - Rsep OMP = {:.2f} mm".format(rminrsep_omp)

# Plotting commands.
fontsize = 16
lw = 5
fig, ax = plt.subplots(figsize=(10,6))
ax.axhline(0, color='k', linestyle='--', lw=3)
ax.plot(s, vz00, color="C0", label="M = 0.0", lw=lw)
ax.plot(s, vi00+vti_final, color="C0", linestyle='--', lw=lw)
ax.plot(s, vz01, color="C1", label="M =-0.1", lw=lw)
ax.plot(s, vi01+vti_final, color="C1", linestyle='--', lw=lw)
ax.plot(s, vz02, color="C2", label="M =-0.2", lw=lw)
ax.plot(s, vi02+vti_final, color="C2", linestyle='--', lw=lw)
ax.plot(s, vz03, color="C3", label="M =-0.3", lw=lw)
ax.plot(s, vi03+vti_final, color="C3", linestyle='--', lw=lw)
ax.plot(s, vz04, color="C4", label="M =-0.4", lw=lw)
ax.plot(s, vi04+vti_final, color="C4", linestyle='--', lw=lw)
ax.set_xlabel("Distance from Inner Target (m)", fontsize=fontsize)
ax.set_ylabel(" Average W Velocity (m/s)", fontsize=fontsize)
ax.tick_params(which='both', labelsize=12)
ax.legend(fontsize=fontsize, loc="upper left")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.text(0.7, 0.1, text_str, fontsize=12, transform=ax.transAxes)
fig.tight_layout()
fig.show()
