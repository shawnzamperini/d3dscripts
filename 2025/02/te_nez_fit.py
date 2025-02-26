import matplotlib.pyplot as plt
import oedge_plots
import numpy as np
from scipy.optimize import curve_fit
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression


# Load OedgePlots of our simulation
ncpath = "/home/zamp/oedge_files/d3d-w-wall-param-scan-v4-ne-1.nc"
op = oedge_plots.OedgePlots(ncpath)

# Load a 1D array of all the ne's and Te's in each cell
nes = op.read_data_2d("KNBS")
tes = op.read_data_2d("KTEBS")

# Load the neon density for each cahrge state
ne_dens = []
for i in range(0, 10):
	ne_dens.append(op.read_data_2d("DDLIMS", charge=i+1))

# Calculate the average charge state in each cell, weighting each by its density
avg_ne_z = np.zeros(len(ne_dens[0]))
weight_sum = np.zeros(len(ne_dens[0]))
ne_dens_sum = np.sum(ne_dens, axis=0)
for i in range(0, len(avg_ne_z)):
	for j in range(0, 10):

		# Z = j+1 (the charge)
		avg_ne_z[i] += (j+1) * ne_dens[j][i]
		weight_sum[i] += ne_dens[j][i]

# Divide each by 10 to finish the average calculation.
avg_ne_z /= weight_sum

# Remove data above Te=50 eV since we're mostly fully ionized by then
mask = tes < 50
nes = nes[mask]
tes = tes[mask]
avg_ne_z = avg_ne_z[mask]

# Remove NaNs
mask = ~np.isnan(avg_ne_z)
nes = nes[mask]
tes = tes[mask]
avg_ne_z = avg_ne_z[mask]

# See if we can derive a scaling law of Z (ne, Te)
def model(x, a, b, c):
	x1, x2 = x
	return x1**a + x2**b + c

xdata = np.array(list(zip(nes, tes)))
popt, pcov = curve_fit(model, xdata.T, avg_ne_z, 
	bounds=((-10, -10, -np.inf), (10, 10, np.inf)))

# Normalizing data
scaler = StandardScaler()
xdata_norm = scaler.fit_transform(xdata)

# Fit and predict model
linreg = LinearRegression().fit(xdata_norm, avg_ne_z)
avg_ne_z_fit = linreg.predict(xdata_norm)


fig, ax = plt.subplots()

ax.scatter(model(xdata.T, *popt), avg_ne_z)
ax.scatter(avg_ne_z_fit, avg_ne_z)
ax.set_xlim([0, 11])
ax.set_ylim([0, 11])

fig.tight_layout()
fig.show()
