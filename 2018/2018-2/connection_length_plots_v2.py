import numpy as np
import get_ts as ts
import openpyxl as xl
import MDSplus as mds
from scipy import interpolate
import matplotlib.pyplot as plt


def returnArray(sheet, lowRange, highRange):
    cells = sheet[lowRange:highRange]
    cells = np.transpose(cells)
    cells = np.reshape(cells, cells.size)
    values = np.array([cell.value for cell in cells])
    return values

conn_wb = xl.load_workbook("Data/conn_lengths.xlsx", data_only=True)
samp_wb = xl.load_workbook("Data/LP_with_fit.xlsx", data_only=True)
conn_sheet = conn_wb.get_sheet_by_name("Sheet1")
samp_sheet = samp_wb.get_sheet_by_name("Fit Te's")

# Get the sampling length data. This is in R-Rsep vs. L.
samp_R        = returnArray(samp_sheet, "M2", "M348") / 100.0
samp_L_A        = returnArray(samp_sheet, "P2", "P348")
samp_L_B        = returnArray(samp_sheet, "Q2", "Q348")
samp_L_C        = returnArray(samp_sheet, "R2", "R348")

# Get the connection length data. This is in R-Rsep_omp vs. L.
conn_RminRsep_omp_A = returnArray(conn_sheet, "A17", "A113")
conn_odf_A          = returnArray(conn_sheet, "B17", "B113")
conn_idf_A          = returnArray(conn_sheet, "C17", "C113")
conn_RminRsep_omp_B = returnArray(conn_sheet, "E17", "E113")
conn_odf_B          = returnArray(conn_sheet, "F17", "F113")
conn_idf_B          = returnArray(conn_sheet, "G17", "G113")
conn_RminRsep_omp_C = returnArray(conn_sheet, "I17", "I113")
conn_odf_C          = returnArray(conn_sheet, "J17", "J113")
conn_idf_C          = returnArray(conn_sheet, "K17", "K113")

# Load the gfile now.
conn = mds.Connection("localhost")
gfile = ts.load_gfile_mds(167196, 3000, connection=conn)

# Create mesh out of the R's and Z's.
Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])

# Get the lcfs data.
Rs_lcfs = gfile["lcfs"][:,0]
Zs_lcfs = gfile["lcfs"][:,1]
# Find where the max and min of Z is and ignore everything to the left of it.
top_index = np.argmax(Zs_lcfs)
bottom_index = np.argmin(Zs_lcfs)
Rs_lcfs = Rs_lcfs[top_index:bottom_index]
Zs_lcfs = Zs_lcfs[top_index:bottom_index]
f_R_lcfs = interpolate.interp1d(Zs_lcfs, Rs_lcfs, assume_sorted=False)


# Get the R and Z of the magnetic axis.
Z_axis = gfile['ZmAxis']
R_axis = gfile['RmAxis']

# Interpolation for psin(R, Z).
Rs_trunc = Rs > R_axis
f_psin = interpolate.Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile["psiRZn"][Rs_trunc])

# Interpolation for R(psin, Z).
f_R = interpolate.Rbf(gfile["psiRZn"][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc])

# Convert samp_R to psin.
Z_loc = -0.188
samp_psin = np.array([])
for R in samp_R:
    tmp_psin = f_psin(R, Z_loc)
    samp_psin = np.append(samp_psin, tmp_psin)

# Then convert each samp_psin to an R at the omp.
samp_R_omp = np.array([])
for psin in samp_psin:
    tmp_R_omp = f_R(psin, Z_axis)
    samp_R_omp = np.append(samp_R_omp, tmp_R_omp)

# Get the R of the separatrix at the omp.
R_omp = f_R_lcfs(Z_axis)

# Get R-Rsep_omp for the probes. Convert from m to cm.
rminrsep_omp = (samp_R_omp - R_omp)

plt.semilogy(rminrsep_omp, samp_L_A, label="A Sampling Length")
plt.semilogy(rminrsep_omp, samp_L_B, label="B Sampling Length")
plt.semilogy(rminrsep_omp, samp_L_C, label="C Sampling Length")
plt.semilogy(conn_RminRsep_omp_B, conn_odf_B, label="Outer Divertor Connection Length")
plt.semilogy(conn_RminRsep_omp_B, conn_idf_B, label="Inner Divertor Connection Length")
plt.legend(prop={"size":18})
plt.xlabel("R-Rsep omp (cm)")
plt.ylabel("Length (cm)")
plt.title("Connection vs. Sampling Lengths")
plt.show()
