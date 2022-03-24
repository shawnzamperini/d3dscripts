import numpy as np
import matplotlib.pyplot as plt


plt.rcParams['font.family'] = 'Century Gothic'

def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


# Assume the average blob velocity
rmrs = np.linspace(0, 10, 50)
v_blob1 = (0, 2500)
v_blob2 = (10, 100)
v_hole1 = (0, -2000)
v_hole2 = (3, 0)

m_blob = (v_blob2[1] - v_blob1[1]) / (v_blob2[0] - v_blob1[0])
m_hole = (v_hole2[1] - v_hole1[1]) / (v_hole2[0] - v_hole1[0])

v_blob = m_blob * (rmrs - v_blob1[0]) + v_blob1[1]
v_hole = m_hole * (rmrs - v_hole1[0]) + v_hole1[1]

v_hole[rmrs > 3] = 0.0

v_blob = np.exp(-rmrs / 5) * 2500 + 300
v_hole = -np.exp(-rmrs / 2) * 3500

fig, ax = plt.subplots()
ax.plot(rmrs, v_blob, label="Blob", color="k", lw=4, linestyle=":")
ax.plot(rmrs, v_hole, label="Hole", color="k", lw=4, linestyle="--")
ax.plot(rmrs, v_blob+v_hole, label="Average", lw=4, color="k", linestyle="-")
ax.set_xlabel("R - Rsep (cm)", fontsize=14)
ax.set_ylabel("Velocity (m/s)", fontsize=14)
ax.legend(fontsize=12)
ax.grid()
ax.tick_params(labelsize=12)
ax.axhline(0.0, color="k", lw=2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()
