from gadata import gadata
import pandas as pd
import numpy as np
import MDSplus
import matplotlib.pyplot as plt


def plot_bolo(shot, tmin=2000, tmax=5000, exclude_chords=[]):
    bolo_dict = {}
    bolo_zdata = None
    conn = MDSplus.Connection('localhost')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for chord in range(1, 25):
        print('Chord {}'.format(chord))
        point_name = 'BOL_U' + format(chord, '02d') + '_P'
        bolo_data = gadata(point_name, shot=shot, connection=conn)
        min_idx = np.where(bolo_data.xdata < tmin)[0][0]
        max_idx = np.where(bolo_data.xdata > tmax)[0][0]
        if bolo_zdata is None:
            bolo_zdata = np.zeros((24, max_idx - min_idx))

        if chord in exclude_chords:
            continue
        #print(bolo_data.zdata[min_idx:max_idx])
        bolo_zdata[chord-1] = bolo_data.zdata[min_idx:max_idx]
    # Flip bolo_zdata upside down so that topmost chord (24) is on the top in the graph.
    bolo_zdata = np.flipud(bolo_zdata)
    ax1.pcolor(bolo_zdata)
    fig.show()

    return bolo_zdata

def plot_avg_bolo(shots, tmin=2000, tmax=5000, exclude_chords=[]):
    conn = MDSplus.Connection('localhost')
    all_bolo = None
    shot_count = 0

    for shot in shots:
        print("Shot {}".format(shot))

        # Get data for this shot.
        bolo_zdata = None
        for chord in range(1, 25):
            point_name = 'BOL_U' + format(chord, '02d') + '_P'
            bolo_data = gadata(point_name, shot=shot, connection=conn)

            # Restrict the data to the desired time range. Only need to find
            # once since each chord shares the same time index.
            if bolo_zdata is None:
                min_idx = np.where(bolo_data.xdata < tmin)[0][0]
                max_idx = np.where(bolo_data.xdata > tmax)[0][0]
                bolo_zdata = np.zeros((24, max_idx - min_idx))

            # Don't include the bad chords.
            if chord in exclude_chords:
                continue

            # Put this chord into the array for this shot.
            bolo_zdata[chord-1] = bolo_data.zdata[min_idx:max_idx]

        # Now flip the bolo_zdata so the topmost chord is the top row (for plotting).
        #bolo_zdata = np.flipud(bolo_zdata)

        # Create the array to hold all the data now.
        if all_bolo is None:
            all_bolo = np.zeros((len(shots), 24, max_idx-min_idx))

        # Put this shot's bolo data into the array.
        all_bolo[shot_count] = bolo_zdata
        shot_count += 1

    # Now that we have all the data for each shot in a 3D array, average along
    # the axis to make this a 2D array of the average data.
    avg_bolo = np.mean(all_bolo, axis=0)
    avg_bolo_x = np.linspace(tmin, tmax, len(avg_bolo[0]))
    avg_bolo_y = np.arange(0, 25, 1)

    # Plot it up.
    fig = plt.figure(figsize=(7,5))
    ax1 = fig.add_subplot(111)
    ax1.pcolor(avg_bolo_x, avg_bolo_y, avg_bolo)
    ax1.set_xlim([tmin, tmax])
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Chord')
    fig.tight_layout()
    fig.show()

    return avg_bolo
