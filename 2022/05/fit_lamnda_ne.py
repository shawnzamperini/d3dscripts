# From the OMFIT netCDF files fit an exponential for lambda_ne.
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


shots = [167196,170843,176508,179900,164262,162766,171432,172409,172413,174212,174305]
shots = [167195]

root = "/Users/zamperini/Documents/d3d_work/files/OMFITprofiles_"

def exp_fit(x, a, b, c):
    return a * np.exp(-b * x) + c

fitdata = {}; plotshots = []
for shot in shots:
    print(shot)
    try:
        nc = netCDF4.Dataset("{}{}_FIT.nc".format(root, shot))
    except:
        print(" Missing NetCDF file.")
        continue
    rho = nc["rho"][:]
    ne = nc["n_e"][:].mean(axis=0)
    sep = np.argmin(np.abs(rho - 1))
    nenorm = ne / ne[sep]

    # Data can get unreliable further out, so just restrict it a bit.
    #end = np.argmin(np.abs(rho - 1.08))
    end = 999

    # Perform fit to only values outside the separatrix.
    popt, pcov = curve_fit(exp_fit, rho[sep:end], nenorm[sep:end], maxfev=5000)
    fitdata[shot] = {"rho":rho[sep:], "nenorm":nenorm[sep:], "parms":popt}
    plotshots.append(shot)


# Plot.
fig, ax1 = plt.subplots()
count = 0
print("Lambda_ne's")
for shot in plotshots:
    fitx = np.linspace(1, 1.2, 20)
    fity = exp_fit(fitx, *fitdata[shot]["parms"])
    ax1.plot(fitdata[shot]["rho"], fitdata[shot]["nenorm"], color="C{}".format(count), label=shot)
    ax1.plot(fitx, fity, color="C{}".format(count), linestyle="--")
    count += 1
    print("{}: {:.3f}".format(shot, 1/fitdata[shot]["parms"][1]))
ax1.legend()
ax1.set_xlabel("Rho")
ax1.set_ylabel("ne_norm")
fig.tight_layout()
fig.show()
