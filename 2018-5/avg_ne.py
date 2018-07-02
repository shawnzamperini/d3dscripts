import gadata
import numpy as np
import MDSplus as mds


shots = [167219, 167237, 167266, 167268, 167535, 167530, 167537, 167353,
         167321, 167408, 167277, 167279, 167481]
#shots = [167268]

conn = mds.Connection('localhost')
lower_time = 2000
upper_time = 4000
ne_dict     = {}
a_dict      = {}
ne_err_dict = {}
a_err_dict  = {}
q95_dict     = {}
q95_err_dict = {}
bt_dict     = {}
bt_err_dict = {}
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

    q95_obj  = gadata.gadata('\Q95', shot, connection=conn)
    times = q95_obj.xdata
    q95s  = q95_obj.zdata
    idx    = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_q95s = np.mean(q95s[idx])
    std_q95s = np.std(q95s[idx])

    bt_obj  = gadata.gadata('\BT', shot, connection=conn)
    times = bt_obj.xdata.data()
    bts  = bt_obj.zdata.data()
    idx    = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_bts = np.mean(bts[idx])
    std_bts = np.std(bts[idx])


    ne_dict[str(shot)]      = avg_nes
    a_dict[str(shot)]       = avg_amins
    ne_err_dict[str(shot)]  = std_nes
    a_err_dict[str(shot)]   = std_amins
    q95_dict[str(shot)]     = avg_q95s
    q95_err_dict[str(shot)] = std_q95s
    bt_dict[str(shot)]      = avg_bts
    bt_err_dict[str(shot)]  = std_bts
