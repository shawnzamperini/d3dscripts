# This script is mainly as a beta for testing the fitting on the L-mode probes.
# In the future it would be ideal to use pandas.

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

def plot_it(x, y, y_err, x_fit, y_fit, b, probe):
    plt.rcParams.update({'font.size': 34})
    plt.errorbar(x=x, y=y, yerr=y_err, fmt="b.", ms=15, capsize=5)
    plt.plot(x_fit, y_fit, "b--", ms=15)
    plt.xlabel(r"$\mathrm{R-R_{sep}\ omp\ (cm)}$")
    plt.ylabel(r"$\mathrm{W\ Areal\ Density\ (10^{15} cm^{-2})}$")
    plt.title(probe + " Fitting")
    plt.ylim([0, np.max(y) * 1.1])
    plt.annotate(r"$\mathrm{\lambda = }$" + "{:.2f} cm".format(1/b),
                 xy=(0.6, 0.6),
                 xycoords='axes fraction',
                 bbox=dict(boxstyle="square", fc="w"))
    plt.show(block=False)

def exp_fit(x, a, b, c):
    return a * np.exp(-b * x) + c

# Get the path to the Excel file.
#print("Select Excel file with probe data...\n")
#tkinter.Tk().withdraw()
#probe_excel_file = askopenfilename()
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

while True:
    probe = input("\nWhich probe would you like to plot (A-ITF, A-OTF, B-ITF, ...)? Enter 'q' to quit.\n  ")
    print()
    probe = probe.upper()
    if probe == 'Q':
        break
    exclude = None
    while True:
        if probe == "Q":
            break
        elif probe == "A-ITF":
            x = r_omp_AITF
            y = w_AITF
            y_err = w_AITF_err
        elif probe == "A-OTF":
            x = r_omp_AOTF
            y = w_AOTF
            y_err = w_AOTF_err
        elif probe == "B-ITF":
            x = r_omp_BITF
            y = w_BITF
            y_err = w_BITF_err
        elif probe == "B-OTF":
            x = r_omp_BOTF
            y = w_BOTF
            y_err = w_BOTF_err
        elif probe == "C-ITF":
            x = r_omp_CITF
            y = w_CITF
            y_err = w_CITF_err
        elif probe == "C-OTF":
            x = r_omp_COTF
            y = w_COTF
            y_err = w_COTF_err

        else:
            print("Incorrect probe entry.\n")
            continue

        guess = (100, 1, 0)
        with warnings.catch_warnings():
            # Turn the RuntimeWarning from curve_fit into an error so we can catch it.
            warnings.simplefilter("error", OptimizeWarning)
            try:
                popt, pcov = curve_fit(exp_fit, x[:exclude], y[:exclude], p0=guess, maxfev=50000)
            except:
                print("Error fitting exponential. Try excluding some points.\n")
            else:
                x_upper = np.max(x)
                x_lower = np.min(x)
                x_fit = np.linspace(x_lower, x_upper, 100)
                y_fit = exp_fit(x_fit, *popt)
                plot_it(x, y, y_err, x_fit, y_fit, popt[1], probe)

        try:
            exclude = -1 * int(input("Exclude the furthest ___ points (enter 'q' to quit/select new probe): "))
        except:
            break
