# coding: utf-8
get_ipython().run_line_magic('run', 'div_core.py')
div['psins']
div['psins'].keys()
div['psins']['avg_psins']
core['psins']['avg_psins']
div['psins']['all_psins']
div['psins'].keys()
div.keys()
div['r']
div['channel']
divr = div['r']['Y']
divz = div['z']['Y']
divz
get_ipython().run_line_magic('run', 'get_ts')
core['psins']['avg_psins']
div['psins']['avg_psins']
div['psins']['avg_psins_err']
div['psins']['avg_omps']
ts_dict = div
import numpy as np
import MDSplus as mds
import numpy as np
import scipy.interpolate as scinter
import math
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
tmin=2500
tmax=4500
tstep=250
tree='EFIT01'
# Get the times at which the Te and ne data was taken for each chord.
times = ts_dict["temp"]["X"]
range_wanted = np.arange(tmin, tmax, tstep)

# Get index where the times start and end for the requested range.
min_index = np.argmax(times>tmin)
max_index = np.argmax(times>tmax)

print("Time range used: [" + str(times[min_index]) + ", " + str(times[max_index]) + "]")

# Get the actual min and max times. Fill array of times for gfiles to request.
ts_min_time = times[min_index]
ts_max_time = times[max_index]
times_for_gfile = np.arange(ts_min_time, ts_max_time, tstep)
num_of_times = times_for_gfile.size

# Parameters needed from ts_dict.
shot = ts_dict["shot"]
ts_Rs = ts_dict["r"]["Y"]
ts_Zs = ts_dict["z"]["Y"]
channels = ts_dict["channel"]["Y"]
temp = ts_dict["temp"]["Y"]
dens = ts_dict["density"]["Y"]

# Create connection here then pass it to gfile function so you don't have to
# create a new connection every time. Saves a lot of time.
conn = mds.Connection("localhost")

# 2D array to hold all the psin values. The index of each row corresponds to
# a certain channel/chord.
all_psins = np.zeros((num_of_times, len(channels)))
all_tes   = np.zeros((num_of_times, len(channels)))
all_nes   = np.zeros((num_of_times, len(channels)))
all_omps  = np.zeros((num_of_times, len(channels)))
row_index = 0

# Get the gfiles one time slice at a time.
for time in times_for_gfile:

    # Load gfile. Cast to int since gfiles don't have decimal times.
    gfile = load_gfile_mds(shot, int(time), tree=tree, connection=conn)

    # Create grid of R's and Z's.
    Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
    Z_axis = gfile['ZmAxis']
    R_axis = gfile['RmAxis']
    Zes = np.copy(gfile['lcfs'][:, 1][13:-17])
    Res = np.copy(gfile['lcfs'][:, 0][13:-17])
    Rs_trunc = Rs > R_axis

    # Interpolation functions of psin(R, Z) and R(psin, Z).
    f_psiN = scinter.Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile['psiRZn'][Rs_trunc])
    f_Romp = scinter.Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
    f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)

    # R value of separatrix at omp.
    rSep_omp = f_Rs(Z_axis)

    # Temporary array to hold psins.
    psins = np.array([])
    tes   = np.array([])
    nes   = np.array([])
    rminrsep_omp = np.array([])
    for index in range(0, len(channels)):

        # The R and Z of each channel.
        tmp_R = ts_Rs[index]
        tmp_Z = ts_Zs[index]

        # This is the psin at this specific time.
        tmp_psin = f_psiN(tmp_R, tmp_Z)

        # R at the omp of this psin.
        tmp_Romp = f_Romp(tmp_psin, Z_axis)

        # This is the corresponding Te at this time.
        te_index = np.argmax(times>time)
        tmp_te   = temp[index][te_index]
        tmp_ne   = dens[index][te_index]

        # Add it to psins/tes, where in order it is the channels.
        psins = np.append(psins, tmp_psin)
        tes = np.append(tes, tmp_te)
        nes = np.append(nes, tmp_ne)

        # The R-Rsep_omp position.
        tmp_romp = tmp_Romp - rSep_omp
        rminrsep_omp = np.append(rminrsep_omp, tmp_romp)

    # Put row into all the psins. So a 2D array, where each row is the psins/Tes
    # of every channel at a time.
    all_psins[row_index] = psins
    all_tes[row_index]   = tes
    all_nes[row_index]   = nes
    all_omps[row_index]  = rminrsep_omp

    row_index = row_index + 1

def organize_2d(all_vals):
    """
    Function to reorganize the 2D array daya of Te, ne and psin into an
    array where each row corresponds to the data from a specific chord
    (normally each row is the data from all the chords at a specific time).
    """
    all_vals_org     = np.zeros((len(channels), num_of_times))
    avg_vals_org     = np.zeros(len(channels))
    avg_vals_org_err = np.zeros(len(channels))
    row_index   = 0
    for chord in range(0, len(channels)):
        chord_vals = np.array([])
        for vals_at_time in all_vals:
            chord_vals_at_time = vals_at_time[chord]
            chord_vals = np.append(chord_vals, chord_vals_at_time)
        avg_chord_vals     = np.mean(chord_vals)
        avg_chord_vals_err = np.std(chord_vals)
        avg_vals_org[row_index]     = avg_chord_vals
        avg_vals_org_err[row_index] = avg_chord_vals_err
        all_vals_org[row_index]     = chord_vals
        row_index = row_index + 1
    return all_vals_org, avg_vals_org, avg_vals_org_err

# Get the organized 2D arrays of the psin, Te and ne data.
all_psins_org, avg_psins_org, avg_psins_org_err = organize_2d(all_psins)
all_tes_org,   avg_tes_org,   avg_tes_org_err   = organize_2d(all_tes)
all_nes_org,   avg_nes_org,   avg_nes_org_err   = organize_2d(all_nes)
all_omps_org,  avg_omps_org,  avg_omps_org_err  = organize_2d(all_omps)

# Put these bad boys into the dictionary under psin to keep it organized.
psins_dict = {"channels":channels,           "all_psins":all_psins_org,
              "avg_psins":avg_psins_org,     "avg_Tes":avg_tes_org,
              "avg_nes":avg_nes_org,         "avg_Tes_err":avg_tes_org_err,
              "avg_nes_err":avg_nes_org_err, "avg_psins_err":avg_psins_org_err,
              "all_Tes":all_tes_org,         "all_nes":all_nes_org,
              "all_omps":all_omps_org,       "avg_omps":avg_omps_org,
              "avg_omps_err":avg_omps_org_err,
              "times":times_for_gfile}

ts_dict["psins"] = psins_dict
return ts_dict
time-3500
time=3500
# Load gfile. Cast to int since gfiles don't have decimal times.
gfile = load_gfile_mds(shot, int(time), tree=tree, connection=conn)

# Create grid of R's and Z's.
Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
Z_axis = gfile['ZmAxis']
R_axis = gfile['RmAxis']
Zes = np.copy(gfile['lcfs'][:, 1][13:-17])
Res = np.copy(gfile['lcfs'][:, 0][13:-17])
Rs_trunc = Rs > R_axis

# Interpolation functions of psin(R, Z) and R(psin, Z).
f_psiN = scinter.Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile['psiRZn'][Rs_trunc])
f_Romp = scinter.Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)
divr
myrs = np.linspace(0, 1.75, 100)
divz
myzs = np.linspace(0, -1.5, 100)
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', '')
fig = plt.figure()
fig.add_subplot(111)
x = myzs
y = f_psiN(1.49, myzs)
myzs
myrs
y = f_psiN(myrs, myzs)
f_psiN(myrs, myzs)
fig.plot(myrs, y)
ax1 = fig.add_subplot(111)
ax1.plot(myrs, y)
ax1.plot(myzs, y)
ax1.clear()
ax1.contourf(myrs, myzs, y)
R, Z = np.meshgrid(myrs, myzs)
y
gfile['psiRZn']
ax1.contourf(gfile['Rs'], gfile['Zs'], gfile['psinRZn'])
ax1.contourf(gfile['R'], gfile['Z'], gfile['psinRZn'])
ax1.contourf(gfile['R'], gfile['Z'], gfile['psiRZn'])
ax1.clear()
ax1.contour(gfile['R'], gfile['Z'], gfile['psiRZn'])
R_axis
Rs
Rs[R_trunc]
Rs[Rs_trunc]
Rs
gfile['R']
gfile['Z']
gfile['psiRZn']
myrs
rs
r
div.keys()
divr
divz
fpsiN(1.49, -1.2419)
f_psiN(1.49, -1.2419)
f_psiN(1.49, -1.2419)?
f_psiN(1.55, -1.2419)?
f_psiN(1.55, -1.2419)
shot = 167196
get_ipython().run_line_magic('run', 'div_core.py')
div['psins']['avg_psins']
pz = -0.18
gfile['psiRZn']
gfile['R']
gfile['R'][25]
gfile['R'][50]
gfile['R'][55]
testr = gfile['R'][55]
testz = gfile['Z'][55]
testpsin = gfile['psiRZn'][55]
fig = plt.figure()
ax1 = fig.plot(gfile['Z'], testpsin, '.')
ax1 = fig.add_subplot(111)
ax1.add_subplot(111)
ax1.plot(gfile['Z'], testpsin, '.')
gfile['psiRZn']
testr = np.full(gfile['Z'].shape, 1.49)
testr
f_psiN(testr, gfile['Z'])
ax1.plot(gfile['Z'], f_psiN(testr, gfile['Z']))
testrs = np.full(gfile['Z'].shape, 1.49)
testr = gfile['R'][55]
testrs = np.full(gfile['Z'].shape, testr)
ax1[1].clear()
del ax1[1]
ax1[1]
ax1.lines[1]
del ax1.lines[1]
ax1.redraw()
ax1.draw()
ax1.clear()
ax1.plot(gfile['Z'], f_psiN(testr, gfile['Z']))
ax1.plot(gfile['Z'], f_psiN(testrs, gfile['Z']))
ax1.plot(gfile['Z'], testpsin, '.')
testrs
testpsin = gfile['psiRZn'][55]
testpsin
testr
testrs
testr
ax1.clear()
R_axis
gfile['psiRZn'][0]
gfile['lcfs']
ax1.plot(gfile['lcfs'])
ax1.clear()
ax1.plot(gfile['lcfs'][0], gfile['lcfs'][1])
ax1.plot(gfile['lcfs'][:,0], gfile['lcfs'][:,1])
fig.show()
plt.show()
ax1.show()
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(gfile['lcfs'][:,0], gfile['lcfs'][:,1])
gfile['psiRZn']
ax1.plot(gfile['R'], gfile['Z'], '.')
ax1.clear()
ax1.plot(gfile['lcfs'][:,0], gfile['lcfs'][:,1])
RR, ZZ = np.meshgrid(gfile['R'], gfile['Z'])
zip(RR, ZZ)
list(zip(RR, ZZ))
points = np.vstack(RR.ravel(), ZZ.ravel())
points = np.vstack([RR.ravel(), ZZ.ravel()])
points
ax1.scatter(RR, ZZ)
gfile['lcfs']
plt.plot(RR, ZZ, gfile['psiRZn'])
plt.plot_surface(RR, ZZ, gfile['psiRZn'])
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(RR, ZZ, gfile['psiRZn'])
x = gfile['R']
y = gfile['Z']
data = gfile['psiRZn']
from scipy.interpolate import RegularGridInterpolator
my_f = RegularGridInterpolator((x, y), data)
my_f
testr
testz
my_f(testr, testz)
my_f((testr, testz))
div
get_ipython().run_line_magic('save', 'during_psi')
get_ipython().run_line_magic('save', 'during_psi 1-9999')
