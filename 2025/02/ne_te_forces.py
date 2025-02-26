import oedge_plots
import matplotlib.pyplot as plt


# Load the 4 cases:
#   1: ZNe = 3, Wall Te = 5 eV
#   2: ZNe = 6, Wall Te = 5 eV
#   3: ZNe = 5, Wall Te = 3 eV
#   4: ZNe = 5, Wall Te = 6 eV
ops = []
for i in range(1, 5):
	ncpath = "/home/zamp/oedge_files/d3d-w-wall-param-scan-v4-test-{}.nc".format(i)
	ops.append(oedge_plots.OedgePlots(ncpath))

# Load the FF and FiG forces. For the FF, we just assume the ion is at rest
# since we are close to the target. Consider different rings.
#   17: Near-SOL
#   65: Lower-upper baffle
ffs = {}
figs = {}
rings = [17, 65] 
w_charge = 5
for i in range(1, 5):
	for ring in rings:
		ffs["{}-".format(i) + str(ring)] = ops[i-1].along_ring(ring, "ff", 
			charge=w_charge, plot_it=False)
		figs["{}-".format(i) + str(ring)] = ops[i-1].along_ring(ring, "fig", 
			charge=w_charge, plot_it=False)
	
fig, axs = plt.subplots(1, 2, figsize=(8, 5))

# Near-SOL
axs[0].plot(ffs["1-{}".format(rings[0])][0], ffs["1-{}".format(rings[0])][1], 
	color="r")
axs[0].plot(ffs["2-{}".format(rings[0])][0], ffs["2-{}".format(rings[0])][1], 
	color="g")
axs[0].plot(ffs["3-{}".format(rings[0])][0], ffs["3-{}".format(rings[0])][1], 
	color="b")
axs[0].plot(ffs["4-{}".format(rings[0])][0], ffs["4-{}".format(rings[0])][1], 
	color="k")
#axs[0].set_xlim([80, 91])
#axs[0].set_ylim([1e-17, 1e-15])
axs[0].set_yscale("symlog", linthresh=1e-17)

# Far-SOL
axs[1].axhline(0.0, color="k")
axs[1].plot(figs["1-{}".format(rings[0])][0], figs["1-{}".format(rings[0])][1]
	+ ffs["1-{}".format(rings[0])][1], color="r")
axs[1].plot(figs["2-{}".format(rings[0])][0], figs["2-{}".format(rings[0])][1]
	+ ffs["2-{}".format(rings[0])][1], color="r", linestyle="--")
axs[1].plot(figs["3-{}".format(rings[0])][0], figs["3-{}".format(rings[0])][1]
	+ ffs["3-{}".format(rings[0])][1], color="b")
axs[1].plot(figs["4-{}".format(rings[0])][0], figs["4-{}".format(rings[0])][1]
	+ ffs["4-{}".format(rings[0])][1], color="b", linestyle="--")
axs[1].set_yscale("symlog", linthresh=1e-18)

fig.tight_layout()
fig.show()
