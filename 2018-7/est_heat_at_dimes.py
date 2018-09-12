import get_ts as ts
import numpy as np
import matplotlib.pyplot as plt


def heat_flux(shot, start_time=2000, end_time=4500):
    ts_dict = ts.run_script(shot, 'divertor', tmin=start_time, tmax=end_time)
    time = ts_dict['temp']['X']
    temp_0 = ts_dict['temp']['Y'][0]   # Right above DiMES.
    temp_7 = ts_dict['temp']['Y'][-1]  # Topmost chord near separatrix (167196, 167536).

    # Get the data in the time range.
    low  = np.where(time > start_time)[0].min()
    high = np.where(time < end_time)[0].max()
    time = time[low:high]
    temp_0 = temp_0[low:high]
    temp_7 = temp_7[low:high]
    avg_temp_0 = np.mean(temp_0)
    avg_temp_7 = np.mean(temp_7)
    print("Chord 0: {:.2f}\nChord 7: {:.2f}".format(avg_temp_0, avg_temp_7))

    dens_0 = ts_dict['density']['Y'][0]   # Right above DiMES.
    dens_7 = ts_dict['density']['Y'][-1]  # Topmost chord near separatrix (167196, 167536).

    dens_0 = dens_0[low:high]
    dens_7 = dens_7[low:high]
    avg_dens_0 = np.mean(dens_0)
    avg_dens_7 = np.mean(dens_7)
    print("Chord 0: {:.2f}\nChord 7: {:.2f}".format(avg_dens_0, avg_dens_7))

    # Plot temp and density at each chord location.
    fig = plt.figure()
    ax1 = fig.add_subplot(211) # Chord 7
    ax1.plot(ts_dict['temp']['X'], ts_dict['temp']['Y'][-1], 'r.')
    ax1.set_xlim([0,6000])
    ax1.set_ylim([0,500])
    ax1.set_ylabel('Te (eV)')
    ax2 = fig.add_subplot(212) # Chord 0
    ax2.plot(ts_dict['temp']['X'], ts_dict['temp']['Y'][0], 'r.')
    ax2.set_xlim([0,6000])
    ax2.set_ylim([0,50])
    ax2.set_ylabel('Te (eV)')
    ax2.set_xlabel('Time (ms)')
    fig.tight_layout()
    fig.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(211) # Chord 7
    ax1.plot(ts_dict['density']['X'], ts_dict['density']['Y'][-1], 'r.')
    ax1.set_xlim([0,6000])
    #ax1.set_ylim([0,500])
    ax1.set_ylabel('ne (m-3)')
    ax2 = fig.add_subplot(212) # Chord 0
    ax2.plot(ts_dict['density']['X'], ts_dict['density']['Y'][0], 'r.')
    ax2.set_xlim([0,6000])
    #ax2.set_ylim([0,50])
    ax2.set_ylabel('ne (m-3)')
    ax2.set_xlabel('Time (ms)')
    fig.tight_layout()
    fig.show()

    return ts_dict
