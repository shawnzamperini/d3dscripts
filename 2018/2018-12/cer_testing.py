from gadata import gadata
import MDSplus as mds
import pretty_plots as pp
import numpy as np


conn = mds.Connection('localhost')

names = ['t10', 't11', 't12', 't13', 't14', 't15', 't49', 't16', 't48', 'v32']
rads  = [222.14, 223.45, 224.78, 226.29, 227.8, 229.35, 229.79, 230.8, 231.13, 231.75]

avg_rotcs1 = np.array([]); shot1 = 167196
avg_rotcs2 = np.array([]); shot2 = 167481
avg_rotcs3 = np.array([]); shot3 = 167353
for name in names:
    tag = 'cerarotc' + name
    print('Loading ' + tag + '...')
    gaobj = gadata(tag, shot1, connection=conn)
    time = gaobj.xdata
    rotc = gaobj.zdata

    idx = np.logical_and(time>2500, time<4500)
    avg_rotc = np.mean(rotc[idx])
    avg_rotcs1 = np.append(avg_rotcs1, avg_rotc)

    gaobj = gadata(tag, shot2, connection=conn)
    time = gaobj.xdata
    rotc = gaobj.zdata

    idx = np.logical_and(time>2500, time<4500)
    avg_rotc = np.mean(rotc[idx])
    avg_rotcs2 = np.append(avg_rotcs2, avg_rotc)

    gaobj = gadata(tag, shot3, connection=conn)
    time = gaobj.xdata
    rotc = gaobj.zdata

    idx = np.logical_and(time>2500, time<4500)
    avg_rotc = np.mean(rotc[idx])
    avg_rotcs3 = np.append(avg_rotcs3, avg_rotc)

fig = pp.pplot(rads, avg_rotcs1, label=str(shot1), fmt='-')
fig = pp.pplot(rads, avg_rotcs3, label=str(shot3), fmt='-', fig=fig, color=10)
fig = pp.pplot(rads, avg_rotcs2, xlabel='Radial Location (cm)', ylabel='Toroidal Rotation', label=str(shot2), fig=fig, color=8, fmt='-')
