import numpy as np
import matplotlib.pyplot as plt
import get_ts as ts
import MDSplus as mds
import openpyxl as xl
from scipy import interpolate


def returnArray(sheet, lowRange, highRange):
    cells = sheet[lowRange:highRange]
    cells = np.transpose(cells)
    cells = np.reshape(cells, cells.size)
    values = np.array([cell.value for cell in cells])
    return values

wb = xl.load_workbook("Data/LP_with_fit.xlsx", data_only=True)
sheet = wb.get_sheet_by_name("LP Data")

R_lp_192p1  = returnArray(sheet, "B3", "B54")
Te_lp_192p1 = returnArray(sheet, "D3", "D54")
R_lp_195p2  = returnArray(sheet, "BG3", "BG59")
Te_lp_195p2 = returnArray(sheet, "BI3", "BI59")

# Convert R_lp to m, and shift inward 16 mm (1.6 cm).
R_lp_192p1 = (R_lp_192p1 - 1.6) / 100.0
R_lp_195p2 = (R_lp_195p2 - 1.6) / 100.0

Z_loc = -0.188
# Do the gfile thing for R to psin.
conn = mds.Connection("localhost")
gfile = ts.load_gfile_mds(167192, 2000, connection=conn)
Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
Z_axis = gfile['ZmAxis']
R_axis = gfile['RmAxis']
Zes = np.copy(gfile['lcfs'][:, 1][13:-17])
Res = np.copy(gfile['lcfs'][:, 0][13:-17])
Rs_trunc = Rs > R_axis
f_Rs = interpolate.interp1d(Zes, Res, assume_sorted=False)
f_psin = interpolate.Rbf(Rs, Zs, gfile["psiRZn"])
f_Romp = interpolate.Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
rSep_omp = f_Rs(Z_axis)

psin_lp_192p1 = np.array([])
omps_lp_192p1 = np.array([])
for R in R_lp_192p1:
    psin = f_psin(R, Z_loc)
    omp  = f_Romp(psin, Z_axis)
    tmp_rmrs_omp = omp - rSep_omp
    psin_lp_192p1 = np.append(psin_lp_192p1, psin)
    omps_lp_192p1 = np.append(omps_lp_192p1, tmp_rmrs_omp)

# Do the gfile thing for R to psin, except for 196 now.
conn = mds.Connection("localhost")
gfile = ts.load_gfile_mds(167196, 2000, connection=conn)
Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
Z_axis = gfile['ZmAxis']
R_axis = gfile['RmAxis']
Zes = np.copy(gfile['lcfs'][:, 1][13:-17])
Res = np.copy(gfile['lcfs'][:, 0][13:-17])
Rs_trunc = Rs > R_axis
f_Rs = interpolate.interp1d(Zes, Res, assume_sorted=False)
f_psin = interpolate.Rbf(Rs, Zs, gfile["psiRZn"])
f_Romp = interpolate.Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
rSep_omp = f_Rs(Z_axis)

psin_lp_195p2 = np.array([])
omps_lp_195p2= np.array([])
for R in R_lp_195p2:
    psin = f_psin(R, Z_loc)
    omp  = f_Romp(psin, Z_axis)
    tmp_rmrs_omp = omp - rSep_omp
    psin_lp_195p2 = np.append(psin_lp_195p2, psin)
    omps_lp_195p2 = np.append(omps_lp_195p2, tmp_rmrs_omp)

# gfile thing again but for 167195.
#gfile = ts.load_gfile_mds(167195, 4440, connection=conn)
#Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
#f_psin = interpolate.Rbf(Rs, Zs, gfile["psiRZn"])
#psin_lp_195p2 = np.array([])
#for R in R_lp_195p2:
#    psin = f_psin(R, Z_loc)
#    psin_lp_195p2 = np.append(psin_lp_195p2, psin)

# Get TS data for 167192
ts192 = ts.run_script(167192, "core", tmin=1800, tmax=4800, tstep=500)

# Get TS data for 167196
ts196 = ts.run_script(167196, "core", tmin=1800, tmax=4800, tstep=500)

# Get the average values.
all_psin_192 = ts192["psins"]["all_psins"]
all_omps_192 = ts192["psins"]["all_omps"]
all_te_192   = ts192["psins"]["all_Tes"]
all_psin_196 = ts196["psins"]["all_psins"]
all_omps_196 = ts196["psins"]["all_omps"]
all_te_196   = ts196["psins"]["all_Tes"]

chords = all_psin_192.shape[0]
x_avg_192 = np.array([])
y_avg_192 = np.array([])
x_avg_196 = np.array([])
y_avg_196 = np.array([])
for chord in range(0, chords):
    #x = all_psin_192[chord]
    x = all_omps_192[chord]
    y = all_te_192[chord]
    x_avg_192 = np.append(x_avg_192, np.mean(x))
    y_avg_192 = np.append(y_avg_192, np.mean(y))
    plt.plot(x*100, y, "rx", alpha = 0.2)
    #plt.plot(x, y, "rx", alpha = 0.2)
    #x = all_psin_196[chord]
    x = all_omps_196[chord]
    y = all_te_196[chord]
    x_avg_196 = np.append(x_avg_196, np.mean(x))
    y_avg_196 = np.append(y_avg_196, np.mean(y))
    plt.plot(x*100, y, "bx", alpha=0.2)
    #plt.plot(x, y, "bx", alpha=0.2)

# Plot the avg TS values.
plt.rcParams.update({'font.size': 32})
plt.plot(x_avg_192*100, y_avg_192, "rx", label="TS 196192 (1800-4800 ms)")
plt.plot(x_avg_196*100, y_avg_196, "bx", label="TS 196196 (1800-4800 ms)")
#plt.plot(x_avg_192, y_avg_192, "rx", label="TS 196192 (1800-4800 ms)")
#plt.plot(x_avg_196, y_avg_196, "bx", label="TS 196196 (1800-4800 ms)")
# Then plot the Lp data on it.
#plt.plot(psin_lp_192p1, Te_lp_192p1, 'r>', label="LP 167192 P1")
#plt.plot(psin_lp_195p2, Te_lp_195p2, 'g>', label="LP 167195 P2")
plt.plot(omps_lp_192p1*100, Te_lp_192p1, 'r>', label="LP 167192 P1")
plt.plot(omps_lp_195p2*100, Te_lp_195p2, 'g>', label="LP 167195 P2")
#plt.xlabel(r"$\mathrm{\phi_n}$")
plt.xlabel(r"$\mathrm{R-R_{sep}\ omp\ (cm)}$")
plt.ylabel(r"$\mathrm{T_e\ (eV)}$")
#plt.axis([0.95, 1.4, 0, 80])
plt.ylim([0, 80])
plt.title("Comparing TS to LP on Similar Shots")
plt.legend(prop={"size":28})
plt.show()
