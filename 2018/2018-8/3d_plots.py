import numpy  as np
import pandas as pd
import matplotlib.pyplot  as plt
from mpl_toolkits.mplot3d import Axes3D


probe = 'D'

# Path to filename with the coordinates of the flux lines. Only D is
# implemented right now.
print("Loading {} file...".format(probe))
if probe == 'D':
    filename = '/mnt/c/Users/Shawn/Dropbox/Connection_Lengths/Raw_Connection_Lengths/'  \
               'struct_3Dwall_noSAS/struct_A_face_D_3Dwall_noSAS.dat'
elif probe == 'U':
    pass

line_coords = pd.read_csv(filename, skiprows=48, sep='\s+', names=['X', 'Y',
                          'Z', 'R', 'Phi'], header=None)

# The Z value signaling the end (D) or beginning (U) of the coordinates of
# each field line.
zmax = -0.18465

# The indices of where each line ends (D) or begins (U).
zmax_idx = line_coords['Z'][line_coords['Z'] == zmax].index.values

# Array to hold Dataframes of each field line.
field_lines = np.zeros(len(zmax_idx), dtype=np.object)

if probe == 'D':

    # Get the values from each index, up until the next one in zmax_idx.
    prev_idx = 0; count = 0
    for next_idx in zmax_idx:
        field_line = line_coords[prev_idx:next_idx]
        field_lines[count] = field_line
        #print("Sorted field line #{}...".format(count+1))
        count += 1
        prev_idx = next_idx

elif probe == 'U':
    pass

# Let's plot 10 evenly spaced lines.
n_lines = 1
plot_idx = [int(x*len(field_lines)/n_lines) for x in range(n_lines)]

plt.style.use('seaborn')
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
for idx in range(len(plot_idx)):
    line = field_lines[idx]
    ax1.plot3D(line['X'], line['Y'], line['Z'], label=str(line['R'].iloc[0]))
    ax1.legend()

# Plot a set of points representing the probe. We know the (R,Z,phi), convert
# them to (x,y,z). z=-0.188, phi = 240
# For D.
if probe == 'D':
    r_coords = np.array([228.1925956, 228.6720378, 229.1515679, 229.6311854,
                         230.1108896, 230.59068,   231.0705561, 231.5505174,
                         232.0305633, 232.5106934, 232.990907,  240])
elif probe == 'U':
    pass

x = 0.01*r_coords*np.cos(215*np.pi/180)
y = 0.01*r_coords*np.sin(215*np.pi/180)
z = np.full(len(x), -0.18465)

ax1.plot3D(x, y, z, lw=5)

#fig.tight_layout()
fig.show()
