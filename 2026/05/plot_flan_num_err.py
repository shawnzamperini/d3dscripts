import matplotlib.pyplot as plt


# Data gathered in 2026/num_drift_scan.xlsx, copied here
flan_drift = [-0.6178, -0.0625, -0.0044, -0.0243, -0.0121, -0.0044, -0.0035, 
	-0.0029, -4.00E-09, -0.0044, -0.0245, -0.2224, -0.6182]
analytic_drift = [-0.621890547, -0.062189055, -0.006218905, -0.024509804,
	-0.012376238, -0.006218905, -0.00498008, -0.003561254, -0.000248756,
	-0.006218905, -0.024875622, -0.223880597, -0.621890547]

fontsize = 16
fig, ax = plt.subplots(figsize=(5,4))
ax.plot([0, -0.65], [0, -0.65], color="k", linestyle="--", zorder=15)
ax.scatter(analytic_drift, flan_drift, s=75, color="tab:red", edgecolors="k", 
	zorder=25)
ax.set_xlabel("Analytic (m/s)", fontsize=fontsize)
ax.set_ylabel("Flan (m/s)", fontsize=fontsize)
ax.grid(alpha=0.3)
ax.text(
    0.45, 0.55, "1:1",
    transform=ax.transAxes,    # axes coordinates (0–1)
    rotation=41,
    ha='center', va='center',
    fontsize=fontsize)
fig.tight_layout()
fig.show()
