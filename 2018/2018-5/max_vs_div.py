import get_lp            as lp
import numpy             as np
import matplotlib.pyplot as plt
from scipy.optimize      import curve_fit


shot   = 167196
m_deut = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
if True:
    lp_dict = lp.get_dict_of_lps(shot)

# Get estimated flux at floor ring using the probes there: 15, 17, 19, 21.
avg_flux_floor_vals = np.array([])
avg_flux_floor_stds = np.array([])
for p in [15, 17, 19, 21]:
    try:
        # Load lp data.
        mylp = lp_dict['probe ' + str(p)]

        # Get indices of time range 2000-4000.
        idx  = np.where(np.logical_and(mylp['time']>=2000, mylp['time']<=4000))

        # Get the density and temp data.
        mydens = np.array(mylp['dens'])[idx]
        mytemp = np.array(mylp['temp'])[idx]

        # Calculate sound speed and then flux.
        mycs   = np.sqrt(2 * mytemp / m_deut)
        myflux = 0.5 * mydens * mycs

        # Calculate average flux at the floor ring.
        tmp_avg_flux     = np.mean(myflux)
        tmp_avg_flux_std = np.std(myflux)
        avg_flux_floor_vals = np.append(avg_flux_floor_vals, tmp_avg_flux)
        avg_flux_floor_stds = np.append(avg_flux_floor_stds, tmp_avg_flux_std)

    except:
        print('No data for probe ' + str(p) + '.')

# Finally the overall average.
avg_flux_floor     = np.mean(avg_flux_floor_vals)
avg_flux_floor_std = np.mean(avg_flux_floor_stds)
print("Average flux at floor: {0:.2e} += {1:.3e}".format(avg_flux_floor, avg_flux_floor_std))


# Get estimated flux at shelf ring exponentially extrapolating the probes inwards
# to the shelf ring location: 23, 25, 27, 29, 31.
avg_flux_shelf_vals = np.array([])
avg_flux_shelf_stds = np.array([])
shelf_rs = np.array([])
avail_ps = 0
for p in [23, 25, 27, 29, 31]:
    try:
        # Load lp data.
        mylp = lp_dict['probe ' + str(p)]
        shelf_rs = np.append(shelf_rs, mylp['rprobe'])

        # Get indices of time range 2000-4000.
        idx  = np.where(np.logical_and(mylp['time']>=2000, mylp['time']<=4000))

        # Get the density and temp data.
        mydens = np.array(mylp['dens'])[idx]
        mytemp = np.array(mylp['temp'])[idx]

        # Calculate sound speed and then flux.
        mycs   = np.sqrt(2 * mytemp / m_deut)
        myflux = 0.5 * mydens * mycs

        # Calculate average flux at the floor ring.
        tmp_avg_flux     = np.mean(myflux)
        tmp_avg_flux_std = np.std(myflux)
        avg_flux_shelf_vals = np.append(avg_flux_shelf_vals, tmp_avg_flux)
        avg_flux_shelf_stds = np.append(avg_flux_shelf_stds, tmp_avg_flux_std)
        avail_ps += 1

    except:
        print('No data for probe ' + str(p) + '.')

# About where the center of the shelf ring is.
shelf_ring_r = 1.425

if avail_ps > 2:
    # Now do the exponential extrapolation.
    def exp_fit(x, a, b):
        return a * np.exp(-b * x)

    avg_flux_shelf_vals *= 10e-17
    popt, pcov = curve_fit(exp_fit, shelf_rs, avg_flux_shelf_vals)
    flux_at_shelf_ring = exp_fit(shelf_ring_r, *popt) * 10e17
elif avail_ps == 2:
    # Or do a linear extrapolation.
    m = (avg_flux_shelf_vals[1] - avg_flux_shelf_vals[0]) / (shelf_rs[1] - shelf_rs[0])
    flux_at_shelf_ring = avg_flux_shelf_vals[0] + m * (shelf_ring_r - shelf_rs[0])
else:
    print("Only " + str(avail_ps) + " shelf probe(s). Cannot extrapolate.")

print("Average flux at shelf: {0:.2e}".format(flux_at_shelf_ring))

def plot_shelf_fit():
    x_fit = np.linspace(shelf_rs.min(), shelf_rs.max(), 100)
    y_fit = exp_fit(x_fit, *popt)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(shelf_rs, avg_flux_shelf_vals, '.')
    ax1.plot(x_fit, y_fit)
    fig.show()
