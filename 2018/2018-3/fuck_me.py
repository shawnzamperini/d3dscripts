import matplotlib.pyplot as plt
from tkinter.filedialog import askopenfilename
import tkinter
import openpyxl as xl
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
import warnings


def returnArray(sheet, lowRange, highRange):
    cells = sheet[lowRange:highRange]
    cells = np.transpose(cells)
    cells = np.reshape(cells, cells.size)
    values = np.array([cell.value for cell in cells])
    return values

def exp_fit(x, a, b, c):
    return a * np.exp(-b * x) + c

probe_excel_file = '/home/shawn/d3dscripts/Data/LModeProbes.xlsx'

# Get the workbook and the sheets.
wb = xl.load_workbook(probe_excel_file, data_only=True)
A2_sheet = wb.get_sheet_by_name("A2")
B2_sheet = wb.get_sheet_by_name("B2")
C2_sheet = wb.get_sheet_by_name("C2")

r_omp_AITF = returnArray(A2_sheet, "A2", "A20")
r_omp_AOTF = returnArray(A2_sheet, "M2", "M20")
r_omp_BITF = returnArray(B2_sheet, "A2", "A22")
r_omp_BOTF = returnArray(B2_sheet, "M2", "M22")
r_omp_CITF = returnArray(C2_sheet, "A2", "A22")
r_omp_COTF = returnArray(C2_sheet, "M2", "M22")
w_AITF = returnArray(A2_sheet, "C2", "C20")
w_AOTF = returnArray(A2_sheet, "I2", "I20")
w_BITF = returnArray(B2_sheet, "C2", "C22")
w_BOTF = returnArray(B2_sheet, "I2", "I22")
w_CITF = returnArray(C2_sheet, "C2", "C22")
w_COTF = returnArray(C2_sheet, "I2", "I22")
w_AITF_err = returnArray(A2_sheet, "J2", "J20")
w_AOTF_err = returnArray(A2_sheet, "B2", "B20")
w_BITF_err = returnArray(B2_sheet, "J2", "J22")
w_BOTF_err = returnArray(B2_sheet, "B2", "B22")
w_CITF_err = returnArray(C2_sheet, "J2", "J22")
w_COTF_err = returnArray(C2_sheet, "B2", "B22")

def get_popt(x, y, exclude):
    guess = (100, 1, 0)
    popt, pcov = curve_fit(exp_fit, x[:exclude], y[:exclude], p0=guess, maxfev=50000)
    return popt

aitf_ex = -4
aotf_ex = -2
bitf_ex = -5
botf_ex = -2
citf_ex = -2
cotf_ex = -2
aitf_popt = get_popt(r_omp_AITF, w_AITF, aitf_ex)
aotf_popt = get_popt(r_omp_AOTF, w_AOTF, aotf_ex)
bitf_popt = get_popt(r_omp_BITF, w_BITF, bitf_ex)
botf_popt = get_popt(r_omp_BOTF, w_BOTF, botf_ex)
citf_popt = get_popt(r_omp_CITF, w_CITF, citf_ex)
cotf_popt = get_popt(r_omp_COTF, w_COTF, cotf_ex)
aitf_lamb = aitf_popt[1] ** -1
aotf_lamb = aotf_popt[1] ** -1
bitf_lamb = bitf_popt[1] ** -1
botf_lamb = botf_popt[1] ** -1
citf_lamb = citf_popt[1] ** -1
cotf_lamb = cotf_popt[1] ** -1

def do_plot(x, y, y_err, popt, c, label):
    # Create exp fit data.
    x_upper = np.max(x)
    x_lower = np.min(x)
    x_fit = np.linspace(x_lower, x_upper, 100)
    y_fit = exp_fit(x_fit, *popt)

    # Plot the data that was fitted.
    plt.rcParams.update({'font.size': 34})
    plt.errorbar(x=x, y=y, yerr=y_err, fmt=c+".", ms=15, capsize=5, label=label)
    plt.plot(x_fit, y_fit, c+"--", ms=15)
    plt.xlabel(r"$\mathrm{R-R_{sep}\ omp\ (cm)}$")
    plt.ylabel(r"$\mathrm{W\ Areal\ Density\ (10^{15} cm^{-2})}$")
    #plt.title(probe + " Fitting")
    plt.ylim([0, np.max(y) * 1.1])

# Plot A's.
do_plot(r_omp_AOTF[:-1], w_AOTF[:-1], w_AOTF_err[:-1], aotf_popt, 'b', 'A-OTF')
plt.annotate(r"$\mathrm{\lambda = }$" + "{:.2f} cm".format(aotf_lamb),
             xy=(0.6, 0.6),
             xycoords='axes fraction',
             bbox=dict(boxstyle="square", fc="blue", color='blue', alpha=0.2))
do_plot(r_omp_AITF, w_AITF, w_AITF_err, aitf_popt, 'r', 'A-ITF')
plt.annotate(r"$\mathrm{\lambda = }$" + "{:.2f} cm".format(aitf_lamb),
             xy=(0.6, 0.4),
             xycoords='axes fraction',
             bbox=dict(boxstyle="square", fc="red", color='red', alpha=0.2))
plt.legend()
plt.show()

# Plot B's.
do_plot(r_omp_BOTF, w_BOTF, w_BOTF_err, botf_popt, 'b', 'B-OTF')
plt.annotate(r"$\mathrm{\lambda = }$" + "{:.2f} cm".format(botf_lamb),
             xy=(0.6, 0.6),
             xycoords='axes fraction',
             bbox=dict(boxstyle="square", fc="blue", color='blue', alpha=0.2))
do_plot(r_omp_BITF, w_BITF, w_BITF_err, bitf_popt, 'r', 'B-ITF')
plt.annotate(r"$\mathrm{\lambda = }$" + "{:.2f} cm".format(bitf_lamb),
             xy=(0.6, 0.4),
             xycoords='axes fraction',
             bbox=dict(boxstyle="square", fc="red", color='red', alpha=0.2))
plt.legend()
plt.show()

# Plot C's.
do_plot(r_omp_COTF, w_COTF, w_COTF_err, cotf_popt, 'b', 'C-OTF')
plt.annotate(r"$\mathrm{\lambda = }$" + "{:.2f} cm".format(cotf_lamb),
             xy=(0.6, 0.6),
             xycoords='axes fraction',
             bbox=dict(boxstyle="square", fc="blue", color='blue', alpha=0.2))
do_plot(r_omp_CITF, w_CITF, w_CITF_err, citf_popt, 'r', 'C-ITF')
plt.annotate(r"$\mathrm{\lambda = }$" + "{:.2f} cm".format(citf_lamb),
             xy=(0.6, 0.4),
             xycoords='axes fraction',
             bbox=dict(boxstyle="square", fc="red", color='red', alpha=0.2))
plt.legend()
plt.show()
