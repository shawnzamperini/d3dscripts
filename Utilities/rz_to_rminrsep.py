import numpy   as np
import MDSplus as mds
from scipy.interpolate import Rbf, interp1d


def load_gfile_mds(shot, time, tree="EFIT04", exact=False,
                   connection=None, tunnel=True, verbal=True):
    """
    This is scavenged from the load_gfile_d3d script on the EFIT repository,
    except updated to run on python3.

    shot:       Shot to get gfile for.
    time:       Time of the shot to load gfile for, in ms.
    tree:       One of the EFIT trees to get the data from.
    exact:      If True will raise error if time does not exactly match any gfile
                times. False will grab the closest time.
    connection: An MDSplus connection to atlas.
    tunnel:     Set to True if accessing outside DIII-D network.

    returns:    The requested gfile as a dictionary.
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
    if verbal:
        print("Loading gfile:")
        print("  Shot: " + str(shot))
        print("  Tree: " + tree)
        print("  Time: " + str(time))
    signal = 'GTIME'
    k = np.argmin(np.abs(connection.get(base + signal).data() - time))
    time0 = int(connection.get(base + signal).data()[k])

    if (time != time0):
        if exact:
            raise RuntimeError(tree + ' does not exactly contain time %.2f'\
                               %time + '  ->  Abort')
        else:
            if verbal:
                print('Warning: Closest time is ' + str(time0) +'.')
                #print('Fetching time slice ' + str(time0))
            time = time0

    # store data in dictionary
    g = {'shot': shot, 'time': time}

    # get header line
    try:
        header = connection.get(base + 'ECASE').data()[k]
    except:
        print("  No header line.")

    # get all signals, use same names as in read_g_file
    translate = {'MW': 'NR', 'MH': 'NZ', 'XDIM': 'Xdim', 'ZDIM': 'Zdim', 'RZERO': 'R0',
                 'RMAXIS': 'RmAxis', 'ZMAXIS': 'ZmAxis', 'SSIMAG': 'psiAxis', 'SSIBRY': 'psiSep',
                 'BCENTR': 'Bt0', 'CPASMA': 'Ip', 'FPOL': 'Fpol', 'PRES': 'Pres',
                 'FFPRIM': 'FFprime', 'PPRIME': 'Pprime', 'PSIRZ': 'psiRZ', 'QPSI': 'qpsi',
                 'NBBBS': 'Nlcfs', 'LIMITR': 'Nwall'}
    for signal in translate:
        try:
            g[translate[signal]] = connection.get(base + signal).data()[k]
        except:
            print("  Node not found: " + base + signal)

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
    for item in ['Xdim', 'Zdim', 'R0', 'R1', 'RmAxis', 'ZmAxis', 'psiAxis',
                 'psiSep', 'Bt0', 'Ip']:
        g[item] = np.round(np.float64(g[item]), 7)

    # convert single arrays (float32) to double arrays (float64)
    for item in ['Fpol', 'Pres', 'FFprime', 'Pprime', 'psiRZ', 'qpsi',
                 'lcfs', 'wall']:
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


def rz_to_rminrsep_omp(r, z, shot, times, conn=None, choose_interp_region=False, tree='EFIT01'):
    """

    """

    # Make sure times is numpy array.
    times = np.array([times])

    # If conn is None just default to a localhost connection.
    if conn == None:
        connection = mds.Connection("localhost")

    for time in times:

        # Load the gfile.
        gfile = load_gfile_mds(shot, time, tree=tree, connection=conn)

        # R's and Z's of the separatrix, in m.
        outboard_lcfs_rs = gfile['lcfs'][:, 0][13:71]
        outboard_lcfs_zs = gfile['lcfs'][:, 1][13:71]

        # R's and Z's of the grid.
        Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
        Z_axis = gfile['ZmAxis']
        R_axis = gfile['RmAxis']

        # Limit ourselves to values on the outboard side of things. Makes fits
        # cooperate.
        Rs_trunc = Rs > R_axis

        # Interpolation functions of psin(R, Z) and R(psin, Z). Rs_trunc
        # helps with not interpolating the entire plasma, and just that
        # to the right of the magnetic axis, which is normally good enough.
        f_psin = Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile['psiRZn'][Rs_trunc])
        f_Romp = Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
        f_Rs = interp1d(outboard_lcfs_rs, outboard_lcfs_zs, assume_sorted=False)

        # Get the psin values of each input r, z.
        psin = f_psin(r, z)

        # Use the psins to find the value at the omp (Z=Z_axis).
        r_omp = f_Romp(psin, np.full(len(psin), Z_axis))

        # Get R of separatrix at omp and calculate R-Rsep OMP.
        rsep_omp = f_Rs(Z_axis)
        rminrsep_omp = r_omp - rsep_omp

        return rminrsep_omp
