import pandas as pd
import numpy  as np
from scipy.interpolate import interpolate
from scipy.optimize    import curve_fit


probe_file = '/home/shawn/d3dscripts/Data/LModeProbes.xlsx'
aprobe_df = pd.read_excel(probe_file, sheet_name = 'A2', usecols="A:N")[:19]
bprobe_df = pd.read_excel(probe_file, sheet_name = 'B2', usecols="A:N")[:21]
cprobe_df = pd.read_excel(probe_file, sheet_name = 'C2', usecols="A:N")[:21]
lp_file = '/home/shawn/d3dscripts/Data/LPData.xlsx'
lp_df = pd.read_excel(lp_file, sheet_name='LP Data', skiprows=[0], usecols='B:G')[:62]


smol_df = pd.DataFrame()
smol_df['R-Rsep U']    = bprobe_df['rminrsep_U']
smol_df['R-Rsep D']    = bprobe_df['rminrsep_D']
smol_df['W Fluence U'] = bprobe_df['w_areal_U']
smol_df['W Fluence D'] = bprobe_df['w_areal_D']
#smol_df['R-Rsep LP']   = lp_df['R-Rsep (cm)']
#smol_df['ne (e18)']    = lp_df['Ne (E18 m-3)']
#smol_df['Te']          = lp_df['Te (eV)']

m_deut = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
lp_df['cs'] = np.sqrt(2 * lp_df['Te (eV)'] / m_deut)

time = 125
f_ne    = interpolate.interp1d(lp_df['R-Rsep (cm)'], lp_df['Ne (E18 m-3)'])
f_cs    = interpolate.interp1d(lp_df['R-Rsep (cm)'], lp_df['cs'])
x2 = lp_df['R-Rsep (cm)']

def get_frac_U(df):
    f_a = interpolate.interp1d(df['rminrsep_U'], df['w_areal_U'])
    x1 = df['rminrsep_U']
    common_x = np.linspace(np.max([x1.min(), x2.min()]), np.min([x1.max(), x2.max()]), 100)
    new_probe = f_a(common_x)
    new_ne    = f_ne(common_x)
    new_cs    = f_cs(common_x)
    frac = 2 * new_probe*1e19 / (new_ne*1e18 * new_cs * time)
    return common_x, frac

def get_frac_D(df):
    f_a = interpolate.interp1d(df['rminrsep_D'], df['w_areal_D'])
    x1 = df['rminrsep_D']
    common_x = np.linspace(np.max([x1.min(), x2.min()]), np.min([x1.max(), x2.max()]), 100)
    new_probe = f_a(common_x)
    new_ne    = f_ne(common_x)
    new_cs    = f_cs(common_x)
    frac = 2 * new_probe*1e19 / (new_ne*1e18 * new_cs * time)
    return common_x, frac

x_au, frac_au = get_frac_U(aprobe_df)
x_bu, frac_bu = get_frac_U(bprobe_df)
x_cu, frac_cu = get_frac_U(cprobe_df)
x_ad, frac_ad = get_frac_D(aprobe_df)
x_bd, frac_bd = get_frac_D(bprobe_df)
x_cd, frac_cd = get_frac_D(cprobe_df)

# Exponential fit to the density.

def exp_fit(x, a, b):
    return a * np.exp(-b * x)

popt, pcov = curve_fit(exp_fit, lp_df['R-Rsep (cm)'][:52] * 1e-2, lp_df['Ne (E18 m-3)'][:52])
popt_te, pcov_te = curve_fit(exp_fit, lp_df['R-Rsep (cm)'] * 1e-2, lp_df['Te (eV)'])
perr = np.sqrt(np.diag(pcov))
perr_te = np.sqrt(np.diag(pcov_te))


te_lcfs = popt_te[0]
lambda_te = 1 / popt_te[1]
lambda_ne = 1 / popt[1]
lambda_gamma = (1 / lambda_ne + 1 / (2*lambda_te))**(-1)
lambda_ne_err = perr[1] / popt[1] * lambda_ne
lambda_te_err = perr_te[1] / popt_te[1] * lambda_te
ne_lcfs = popt[0]*1e18
dperp = 1
avg_frac = 1e-6
cs_lcfs = np.sqrt(2 * te_lcfs / m_deut)
other_perp_flux_lcfs = 0.5 * avg_frac * cs_lcfs * ne_lcfs

perp_flux_lcfs = avg_frac * dperp * ne_lcfs / lambda_ne
r_range = np.linspace(0, 15) * 1e-2
perp_flux = perp_flux_lcfs * np.exp(-r_range / lambda_gamma)
other_perp_flux = other_perp_flux_lcfs * np.exp(-r_range / lambda_gamma)
other_perp_flux_min = other_perp_flux_lcfs * np.exp(-r_range / (lambda_ne-lambda_ne_err))
other_perp_flux_max = other_perp_flux_lcfs * np.exp(-r_range / (lambda_ne+lambda_ne_err))
