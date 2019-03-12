import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


filename = "track0.txt"

with open(filename) as f:
    contents = f.read()

# Get into entries of strings containing the three data points.
str_entries = contents.split("Particle #7")[1].split('Particle #8')[0].split('\n')[1:-3]

points = np.zeros((len(str_entries), 3))
for i in range(0, len(str_entries)):
    entry = str_entries[i]
    arr = np.array(entry.split(','), dtype=np.float64)
    points[i] = arr

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(points[:,2], points[:,1], -points[:,0], '.')
ax.set_zlim([-0.09, 0])
ax.set_ylim([-0.20, 0.20])
ax.set_xlim([-10, 10])
ax.set_xlabel('Parallel (m)')
ax.set_ylabel('Poloidal (m)')
ax.set_zlabel('Radial (m)')
fig.show()
