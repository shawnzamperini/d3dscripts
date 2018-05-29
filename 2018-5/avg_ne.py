import gadata
import numpy as np
import MDSplus as mds


shots = [167196, 167536, 167534, 167405, 167530, 167463, 167481, 167321,
         167279, 167353, 167277, 167322, 167320]

conn = mds.Connection('localhost')
lower_time = 2000
upper_time = 4000
ne_dict     = {}
a_dict      = {}
ne_err_dict = {}
a_err_dict  = {}
for shot in shots:
    print("Getting values for shot: {}".format(shot))
    lower_time = 2000
    upper_time = 4000
    if shot == 167463:
        upper_time = 3800
    elif shot == 167279:
        lower_time = 3000
        upper_time = 4500
    ne_obj = gadata.gadata('\DENSITY', shot, connection=conn)
    times  = ne_obj.xdata
    nes    = ne_obj.zdata
    idx    = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_nes = np.mean(nes[idx])
    std_nes = np.std(nes[idx])

    a_obj  = gadata.gadata('\AMINOR', shot, connection=conn)
    times = a_obj.xdata
    amins = a_obj.zdata
    idx    = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_amins = np.mean(amins[idx])
    std_amins = np.std(amins[idx])


    ne_dict[str(shot)] = avg_nes
    a_dict[str(shot)]  = avg_amins
    ne_err_dict[str(shot)] = std_nes
    a_err_dict[str(shot)]  = std_amins
