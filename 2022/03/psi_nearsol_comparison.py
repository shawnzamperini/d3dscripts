# Script to compare near-SOL profile of C13 density from our DIVIMP runs.
import oedge_plots
import matplotlib.pyplot as plt


unf_path = "/Users/zamperini/Documents/d3d_work/divimp_files/184527/d3d-184527-inj-013.nc"
fav_path = "/Users/zamperini/Documents/d3d_work/divimp_files/184267/d3d-184267-inj-015.nc"
unf = oedge_plots.OedgePlots(unf_path)
fav = oedge_plots.OedgePlots(fav_path)

# Along ring data. 13 = first SOL ring.
# 184267: 20 has peak Te values at target.
# 184527: Te peaks right from the start at ring 13.
near_ring = 16
sunfnear, nzunfnear = unf.along_ring(near_ring, "DDLIMS", charge="all", plot_it=False)
sfavnear, nzfavnear = fav.along_ring(near_ring, "DDLIMS", charge="all", plot_it=False)
far_ring = 36
sunffar, nzunffar = unf.along_ring(far_ring, "DDLIMS", charge="all", plot_it=False)
sfavfar, nzfavfar = fav.along_ring(far_ring, "DDLIMS", charge="all", plot_it=False)

# Distance from separatrix of this ring.
rminrsep_omp = unf.nc.variables["MIDIST"][1][near_ring]
print("Near: R-Rsep OMP = {:.2f} mm".format(rminrsep_omp * 1000))
rminrsep_omp = unf.nc.variables["MIDIST"][1][far_ring]
print("Far: R-Rsep OMP = {:.2f} mm".format(rminrsep_omp * 1000))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

ax1.plot(sunfnear, nzunfnear, color="tab:purple", label="Unfavorable")
ax1.plot(sfavnear, nzfavnear, color="tab:red", label="Favorable")
ax1.set_xlabel("Distance from outer target (m)")
ax1.set_ylabel("C13 Density (m-3)")
ax1.legend()
ax1.set_ylim(0, 8e15)

ax2.plot(sunffar, nzunffar, color="tab:purple", label="Unfavorable")
ax2.plot(sfavfar, nzfavfar, color="tab:red", label="Favorable")
ax2.set_xlabel("Distance from outer target (m)")
ax2.set_ylabel("C13 Density (m-3)")
ax2.legend()
ax2.set_ylim(0, 8e15)

fig.tight_layout()
fig.show()
