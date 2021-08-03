import oedge_plots
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Some constants.
#ring = 40
mz_amu = 183.84
md_amu = 2.01
mz = mz_amu * 1.66e-27  # Mass of W in kg
md = md_amu * 1.66e-27   # Mass of D in kg
col_log = 15

unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031a.nc"
fav_flow_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006d.nc"
fav_noflow_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc"
unf = oedge_plots.OedgePlots(unf_path)
fav_noflow = oedge_plots.OedgePlots(fav_noflow_path)
fav_flow = oedge_plots.OedgePlots(fav_flow_path)

fav_ring = 18
unf_ring = 70

s_unf, ti_unf = unf.along_ring(unf_ring, "KTIBS", plot_it=False)
s_fav, ti_fav = fav_flow.along_ring(fav_ring, "KTIBS", plot_it=False)
s_unf, ne_unf = unf.along_ring(unf_ring, "KNBS", plot_it=False)
s_fav, ne_fav = fav_flow.along_ring(fav_ring, "KNBS", plot_it=False)
s_fav_no, ti_fav_no = fav_noflow.along_ring(fav_ring, "KTIBS", plot_it=False)
s_fav_no, ne_fav_no = fav_noflow.along_ring(fav_ring, "KNBS", plot_it=False)
#charges = np.arange(1, 30+1)
charges = [8]
#ti_unf = ti_unf * 2
#ti_fav = ti_fav * 2

# The parallel background velocities (vi) for each case.
s_unf, vi_unf  = unf.along_ring(unf_ring, 'Velocity', plot_it=False)
s_fav, vi_fav  = fav_flow.along_ring(fav_ring, 'Velocity', plot_it=False)
s_fav_no, vi_fav_no  = fav_noflow.along_ring(fav_ring, 'Velocity', plot_it=False)

unf_velavg_finals = np.zeros(len(s_unf))
fav_velavg_finals = np.zeros(len(s_fav))
fav_no_velavg_finals = np.zeros(len(s_fav_no))

# Loop for each op case.
count = 0
for op in [unf, fav_flow, fav_noflow]:

    # Choose correct ring.
    if op == unf:
        print("Unf ring")
        ring = unf_ring
        s = s_unf
        ti = ti_unf
        ne = ne_unf
    else:
        print("Fav ring")
        ring = fav_ring
        s = s_fav
        ti = ti_fav
        ne = ne_fav

    # Arrays to hold results for each charge state.
    vti_all    = np.zeros((len(charges), len(s)))
    velavg_all = np.zeros((len(charges), len(s)))
    gamma_all  = np.zeros((len(charges), len(s)))
    nz_all     = np.zeros((len(charges), len(s)))

    for i in range(0, len(charges)):
        charge = charges[i]

        # Get forces.
        #s, ff = op.along_ring(ring, 'ff', plot_it=False, charge=charge)
        s, fig = op.along_ring(ring, 'fig', plot_it=False, charge=charge)

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

    # Average of velavg weighted by densities of each charge state.
    #keep = velavg_all != 0
    #velavg_final = np.average(velavg_all[keep], axis=0, weights=nz_all[keep].mean(axis=0))

    # Store for this case.
    if op == unf:
        unf_velavg_finals = velavg_final
        unf_gamma_finals = gamma_final
        unf_nz_finals = nz_final
        unf_vti_finals = vti_all.mean(axis=0)
    elif op == fav:
        fav_velavg_finals = velavg_final
        fav_gamma_finals = gamma_final
        fav_nz_finals = nz_final
        fav_vti_finals = vti_all.mean(axis=0)
    elif op == fav_noflow:
        fav_no_velavg_finals = velavg_final
        fav_no_gamma_finals = gamma_final
        fav_no_nz_finals = nz_final
        fav_no_vti_finals = vti_all.mean(axis=0)
    #velavg_finals[count] = velavg_final
    count += 1

# This won't change for any of the cases so do it after everyhting.
#vti_final = vti_all.mean(axis=0)

# Index just for shorter variable names. Fix when more cases are added.
#vz_unf = velavg_finals[0]
#vz_fav = velavg_finals[1]

# Get R-Rsep OMP.
rminrsep_omp = unf.nc.variables["MIDIST"][:].data[1][unf_ring] * 1000
text_str = "R - Rsep OMP = {:.2f} mm".format(rminrsep_omp)
print("Unf: {}".format(text_str))
rminrsep_omp = fav.nc.variables["MIDIST"][:].data[1][fav_ring] * 1000
text_str = "R - Rsep OMP = {:.2f} mm".format(rminrsep_omp)
print("Fav: {}".format(text_str))

# Colors from magma colormap.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))

fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.axhline(0.0, linestyle="-", color="k")
#ax.plot(s_unf, unf_velavg_finals)
#ax.plot(s_unf, vi_unf + unf_vti_final, linestyle="--")
#ax.plot(s_unf, vi_unf, linestyle="--")
#ax.plot(s_unf, unf_vti_final, linestyle="--")
ax1.plot(s_fav, fav_velavg_finals, label="vz")
ax1.plot(s_fav, vi_fav + fav_vti_finals, linestyle="--", label="vi+vTi")
ax1.plot(s_fav, vi_fav, linestyle="--", label="vi")
ax1.plot(s_fav, fav_vti_finals, linestyle="--", label="vTi")
ax1.legend()

ax2.plot(s_fav, fav_velavg_finals / np.nanmax(fav_velavg_finals), label="vz")
ax2.plot(s_fav, fav_gamma_finals / np.nanmax(fav_gamma_finals), label="gamma")
ax2.plot(s_fav, fav_nz_finals / np.nanmax(fav_nz_finals), label="nz")
ax2.legend()

fig.tight_layout()
fig.show()
