import matplotlib.pyplot as plt
import numpy as np
import pickle


shots = [200083, 201661, 203188, 200682, 193765, 170361]

times = {200083: 2400, 201661: 3000, 203188: 3000, 200682: 3500, 
	193765: 5500, 170361: 4300}
colors = {200083: "tab:blue", 201661: "tab:red", 203188: "tab:orange", 200682: "tab:pink", 
	193765: "tab:cyan", 170361: "tab:brown"}

root = "/home/zamp/data/"

fig, ax1 = plt.subplots(figsize=(6, 9))
ax1.set_aspect("equal")

for shot in shots:

	with open(root + str(shot) + "_" + str(times[shot]) + ".pickle", "rb") as f:
		gfile = pickle.load(f)	
	
	# Plot separatrices
	ax1.plot(gfile["RBBBS"], gfile["ZBBBS"], label=shot, lw=2, 
		color=colors[shot])
	R, Z = np.meshgrid(gfile["R"], gfile["Z"])
	ax1.contour(R, Z, gfile["PSIRZ_NORM"], levels=[1.0], 
		linewidths=3.0, colors=colors[shot])

# Plot the NT wall
with open(root + str(193765) + "_" + str(times[193765]) + ".pickle", "rb") as f:
	gfile = pickle.load(f)	
	ax1.plot(gfile["RLIM"], gfile["ZLIM"], color=colors[193765], lw=3, linestyle="--")

# Plot the most recent wall
with open(root + str(203188) + "_" + str(times[203188]) + ".pickle", "rb") as f:
	gfile = pickle.load(f)	
	ax1.plot(gfile["RLIM"], gfile["ZLIM"], color="k", lw=5)

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
leg = ax1.legend(loc="upper center", fontsize=14, ncols=2)
leg.get_frame().set_alpha(1.0)
ax1.set_ylim([-2, 2])

fig.tight_layout()
fig.show()
