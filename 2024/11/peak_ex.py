import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks

# Generate sample data (2D Gaussian peaks)
x = np.linspace(-10, 10, 100)
y = np.linspace(-10, 10, 100)
x, y = np.meshgrid(x, y)
z = np.exp(-(x**2 + y**2)) + np.exp(-((x-5)**2 + (y-5)**2))

# Apply Gaussian filter to smooth the data
z_smooth = gaussian_filter(z, sigma=1)

# Find peaks in the smoothed data
peaks = []
for i in range(z_smooth.shape[0]):
    row_peaks, _ = find_peaks(z_smooth[i, :])
    for peak in row_peaks:
        peaks.append((i, peak))

# Plot the data and detected peaks
plt.imshow(z_smooth, extent=(-10, 10, -10, 10), origin='lower')
plt.colorbar()
for peak in peaks:
    plt.plot(peak[1] - 10, peak[0] - 10, 'rx')  # Adjust coordinates for plotting
plt.title('2D Peak Detection')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()

