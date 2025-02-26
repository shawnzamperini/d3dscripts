import numpy as np
import matplotlib.pyplot as plt


# Root directory. Must mount Drive with command:
# sudo mount -t drvfs G: /mnt/g
root = "/mnt/g/My Drive/Research/Data/sicf_sic_heat_flux/"

# Shots the cap was in for
shots = np.arange(200684, 200689, dtype=int)

# Load deltaR and rdat
r = np.loadtxt("{}/rdat.txt".format(root))
dr = np.loadtxt("{}/deltaR_prof.txt".format(root))

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

for shot in shots:
	
	# Load mean and std heat flux data - plotte against deltaR
	avg_path = "{}/q_prof_mean_{}.txt".format(root, shot)
	std_path = "{}/q_prof_std_{}.txt".format(root, shot)
	q_avg = np.loadtxt(avg_path)
	q_std = np.loadtxt(std_path)

	# Load the inter-ELM heat flux - plotted against rdat
	avg_noelm_path = "{}/q_r_avg_noELMs_{}.txt".format(root, shot)
	q_avg_noelm = np.loadtxt(avg_noelm_path)

	# Add to plot
	axs[0].fill_between(dr, q_avg-q_std, q_avg+q_std, alpha=0.25)
	axs[0].plot(dr, q_avg, label=shot)
	axs[1].plot(r, q_avg_noelm, label=shot)

axs[0].legend()
fig.tight_layout()
fig.show()
