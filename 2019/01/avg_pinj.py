from gadata import gadata
import MDSplus as mds
import numpy as np
import pretty_plots as pp


def run(shot):
    conn = mds.Connection('localhost')
    ga_obj = gadata('pinj', shot, connection=conn)
    time = ga_obj.xdata
    pinj = ga_obj.zdata

    fig = pp.pplot(time, pinj, fmt='-', xlabel='Time (ms)', ylabel='PINJ (kW)', xrange=[0,6000])
    minmax = input('Enter time min/max for analysis range (separated by commas): ').split(',')

    min_idx = np.where(time > float(minmax[0]))[0][0]
    max_idx = np.where(time > float(minmax[1]))[0][0]

    avg_pinj = np.mean(pinj[min_idx:max_idx])

    # Convert from kW to MW.
    print('Average PINJ: {:.2f} MW'.format(avg_pinj/1000))
