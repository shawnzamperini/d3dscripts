import MDSplus as mds
from MDSplus import MdsIpException
import numpy as np
from scipy import interpolate

def load_gfile_mds(shot, time, tree="EFIT01", exact=False, connection=None, tunnel=True):
    """
    This is scavenged from the load_gfile_d3d script on the EFIT repository,
    except updated it to run on python3. Also took out the write2file part since
    I never use it.
    """

    # Connect to server, open tree and go to g-file
    if connection is None:
        if tunnel is True:
            connection = mds.Connection("localhost")
        else:
            connection = mds.Connection('atlas.gat.com')
    connection.openTree(tree, shot)

    base = 'RESULTS:GEQDSK:'

    # get time slice
    print("\nLoading gfile:")
    print("  Shot: " + str(shot))
    print("  Tree: " + tree)
    print("  Time: " + str(time))
    signal = 'GTIME'
    k = np.argmin(np.abs(connection.get(base + signal).data() - time))
    time0 = int(connection.get(base + signal).data()[k])

    if (time != time0):
        if exact:
            raise RuntimeError(tree + ' does not exactly contain time %.2f' %time + '  ->  Abort')
        else:
            print('Warning: ' + tree + ' does not exactly contain time %.2f' %time + ' the closest time is ' + str(time0))
            print('Fetching time slice ' + str(time0))
            time = time0

    # store data in dictionary
    g = {'shot': shot, 'time': time}

    # get header line
    header = connection.get(base + 'ECASE').data()[k]

    # get all signals, use same names as in read_g_file
    translate = {'MW': 'NR', 'MH': 'NZ', 'XDIM': 'Xdim', 'ZDIM': 'Zdim', 'RZERO': 'R0',
                 'RMAXIS': 'RmAxis', 'ZMAXIS': 'ZmAxis', 'SSIMAG': 'psiAxis', 'SSIBRY': 'psiSep',
                 'BCENTR': 'Bt0', 'CPASMA': 'Ip', 'FPOL': 'Fpol', 'PRES': 'Pres',
                 'FFPRIM': 'FFprime', 'PPRIME': 'Pprime', 'PSIRZ': 'psiRZ', 'QPSI': 'qpsi',
                 'NBBBS': 'Nlcfs', 'LIMITR': 'Nwall'}
    for signal in translate:
        g[translate[signal]] = connection.get(base + signal).data()[k]

    g['R1'] = connection.get(base + 'RGRID').data()[0]
    g['Zmid'] = 0.0

    RLIM = connection.get(base + 'LIM').data()[:, 0]
    ZLIM = connection.get(base + 'LIM').data()[:, 1]
    g['wall'] = np.vstack((RLIM, ZLIM)).T

    RBBBS = connection.get(base + 'RBBBS').data()[k][:int(g['Nlcfs'])]
    ZBBBS = connection.get(base + 'ZBBBS').data()[k][:int(g['Nlcfs'])]
    g['lcfs'] = np.vstack((RBBBS, ZBBBS)).T

    KVTOR = 0
    RVTOR = 1.7
    NMASS = 0
    RHOVN = connection.get(base + 'RHOVN').data()[k]

    # convert floats to integers
    for item in ['NR', 'NZ', 'Nlcfs', 'Nwall']:
        g[item] = int(g[item])

    # convert single (float32) to double (float64) and round
    for item in ['Xdim', 'Zdim', 'R0', 'R1', 'RmAxis', 'ZmAxis', 'psiAxis', 'psiSep', 'Bt0', 'Ip']:
        g[item] = np.round(np.float64(g[item]), 7)

    # convert single arrays (float32) to double arrays (float64)
    for item in ['Fpol', 'Pres', 'FFprime', 'Pprime', 'psiRZ', 'qpsi', 'lcfs', 'wall']:
        g[item] = np.array(g[item], dtype=np.float64)

    # Construct (R,Z) grid for psiRZ
    g['dR'] = g['Xdim']/(g['NR'] - 1)
    g['R'] = g['R1'] + np.arange(g['NR'])*g['dR']

    g['dZ'] = g['Zdim']/(g['NZ'] - 1)
    NZ2 = int(np.floor(0.5*g['NZ']))
    g['Z'] = g['Zmid'] + np.arange(-NZ2, NZ2+1)*g['dZ']

    # normalize psiRZ
    g['psiRZn'] = (g['psiRZ'] - g['psiAxis']) / (g['psiSep'] - g['psiAxis'])

    return g

def return_avg(shots, start_time=2500, end_time=5000, time_step=500, tree='EFIT01', server='localhost'):
    z_A = -0.18
    z_B = -0.1546
    z_C = -0.2054
    print()
    print("Z Location of A Probe: " + str(z_A) + " m")
    print("Z Location of B Probe: " + str(z_B) + " m")
    print("Z Location of C Probe: " + str(z_C) + " m")

    # Array to hold al lthe values.
    z_arr      = np.array([])
    rsep_arr_A = np.array([])
    rsep_arr_B = np.array([])
    rsep_arr_C = np.array([])

    # MDSplus connection.
    try:
        conn = mds.Connection(server)

        # Make sure shots was entered in list format. If it's a single int/float,
        # then just put it in a list by itself.
        if isinstance(shots, int) or isinstance(shots, float):
            shots = [shots]

            # Go through one shot at a time.
            for shot in shots:
                time = start_time
                while time <= end_time:
                    # Load the gfile for shot, time.
                    gfile = load_gfile_mds(shot, time, tree='EFIT01', connection=conn)

                    # Z of magnetic axis.
                    z_axis = gfile['ZmAxis']

                    # Z's of the separatrix.
                    Zes = np.copy(gfile['lcfs'][:, 1][13:-12])
                    # R's of the separatrix.
                    Res = np.copy(gfile['lcfs'][:, 0][13:-12])
                    # Interpolation function to get Rsep when given a Z.
                    f_Rs = interpolate.interp1d(Zes, Res, assume_sorted=False)

                    # Rsep for each z location of the three probes.
                    r_sep = {}
                    r_sep['a'] = f_Rs(-0.18)
                    r_sep['b'] = f_Rs(-0.1546)
                    r_sep['c'] = f_Rs(-0.2054)

                    # Put the values into our arrays.
                    z_arr = np.append(z_arr, z_axis)
                    rsep_arr_A = np.append(rsep_arr_A, r_sep['a'])
                    rsep_arr_B = np.append(rsep_arr_B, r_sep['b'])
                    rsep_arr_C = np.append(rsep_arr_C, r_sep['c'])

                    # Go to next time and get the next gfile.
                    time += time_step

            # Calculate average and std. dev. of Z of magnetic axis.
            avg_z          = np.mean(z_arr)
            std_dev_z      = np.std(z_arr)
            num_of_samples = len(shots) * int((end_time-start_time)/time_step)

            # Is error divided by sqrt(N), or just std. dev.?
            #err_z          = std_dev_z / math.sqrt(num_of_samples)
            err_z          = std_dev_z

            # Likewise for the Rsep at the probes.
            avg_rsep_A = np.mean(rsep_arr_A)
            avg_rsep_B = np.mean(rsep_arr_B)
            avg_rsep_C = np.mean(rsep_arr_C)
            #avg_rsep_A_err = np.std(rsep_arr_A) / np.sqrt(num_of_samples)
            #avg_rsep_B_err = np.std(rsep_arr_B) / np.sqrt(num_of_samples)
            #avg_rsep_C_err = np.std(rsep_arr_C) / np.sqrt(num_of_samples)
            avg_rsep_A_err = np.std(rsep_arr_A)
            avg_rsep_B_err = np.std(rsep_arr_B)
            avg_rsep_C_err = np.std(rsep_arr_C)

            returned_info = {"Average Z Magnetic Axis":avg_z,      "Z Error":err_z,
                             "Average Rsep at A Probe":avg_rsep_A, "A Error":avg_rsep_A_err,
                             "Average Rsep at B Probe":avg_rsep_B, "B Error":avg_rsep_B_err,
                             "Average Rsep at C Probe":avg_rsep_C, "C Error":avg_rsep_C_err}

            return returned_info

        else:
            print("Incorrect shot entry")

    except MdsIpException:
        print("Please ssh into atlas.gat.com first.")
