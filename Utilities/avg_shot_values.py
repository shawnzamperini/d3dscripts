from gadata import gadata
import MDSplus as mds
import numpy as np
import pretty_plots as pp


def run(shot, tag, units, mult=1):
    conn = mds.Connection('localhost')
    ga_obj = gadata(tag, shot, connection=conn)
    time = ga_obj.xdata
    value = ga_obj.zdata * mult

    fig = pp.pplot(time, value, fmt='-', xlabel='Time (ms)', ylabel=tag.upper() + ' (' + units + ')', xrange=[0,6000])
    minmax = input('Enter time min/max for analysis range (separated by commas): ').split(',')

    min_idx = np.where(time > float(minmax[0]))[0][0]
    max_idx = np.where(time > float(minmax[1]))[0][0]

    avg_value = np.mean(value[min_idx:max_idx])

    # Convert from kW to MW.
    print('Average {}: {:.3f} {}'.format(tag.upper(), avg_value, units))

    return ga_obj

def ip(shot):
    ga_obj = run(shot, 'IP', 'MA', mult=1e-6)
    return ga_obj

def bt(shot):
    ga_obj = run(shot, 'BT', 'T')
    return ga_obj

def rvsout(shot):
    ga_obj = run(shot, 'RVSOUT', 'm')
    return ga_obj

def R0(shot):
    ga_obj = run(shot, 'R0', 'm')
    return ga_obj

def a(shot):
    ga_obj = run(shot, 'AMINOR', 'm')
    return ga_obj

def kappa(shot):
    ga_obj = run(shot, 'KAPPA', '')
    return ga_obj

def density(shot):
    ga_obj = run(shot, 'DENSITY', 'e13 m-3', mult=1e-13)
    return ga_obj

def pinj(shot):
    ga_obj = run(shot, 'PINJ', 'kW')
    return ga_obj

def prad(shot):
    ga_obj = run(shot, 'PRAD_TOT', 'W')
    return ga_obj
