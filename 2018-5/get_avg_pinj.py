import MDSplus as mds
import numpy   as np


conn = mds.Connection('localhost')

shots = [167196, 167536, 167534, 167405, 167530, 167463, 167481, 167321,
         167279, 167353, 167277, 167322, 167320]

pinj_dict  = {}
ip_dict    = {}
nebar_dict = {}
ip_err_dict    = {}
nebar_err_dict = {}
path_pinj  = '\D3D::TOP.NB:PINJ'
path_ip    = '\D3D::TOP.MHD.IP_PROBES:IP_PROBES'
path_nebar = '\D3D::TOP.ELECTRONS.CO2:NELINE_V2'
for shot in shots:
    print("Calculating average values for shot: " + str(shot))

    # Open the d3d tree.
    tree  = conn.openTree('d3d', shot)

    # Time ranges to get average values for.
    lower_time = 2000
    upper_time = 4000

    # Density spike at the end. Not sure what to make of it.
    if shot == 167463:
        upper_time = 3800

    # Calculate the average PINJs.
    pinjs = conn.get(path_pinj).data()
    times = conn.get('DIM_OF('+path_pinj+')').data()
    idx = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_pinj = np.mean(pinjs[idx])
    std_pinj = np.std(pinjs[idx])

    # Calculate the average IPs.
    ips      = conn.get(path_ip).data()
    times_ip = conn.get('DIM_OF('+path_ip+')').data()
    idx = np.where(np.logical_and(times_ip>=lower_time, times_ip<=upper_time))
    avg_ip   = np.mean(ips[idx])
    std_ip   = np.std(ips[idx])

    # Calculate ne_bar (line integrated electron density).
    nes       = conn.get(path_nebar).data()
    times_nes = conn.get('DIM_OF('+path_nebar+')').data()
    idx = np.where(np.logical_and(times_nes>=lower_time, times_nes<=upper_time))
    avg_nes   = np.mean(nes[idx])
    std_nes   = np.std(nes[idx])

    # Put into dictionaries for use in command line.
    pinj_dict[str(shot)]  = avg_pinj
    ip_dict[str(shot)]    = avg_ip
    nebar_dict[str(shot)] = avg_nes
    ip_err_dict[str(shot)]    = std_ip
    nebar_err_dict[str(shot)] = std_nes
