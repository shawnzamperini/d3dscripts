import pandas as pd
import numpy  as np
import matplotlib.pyplot as plt
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

# The later work. More correct here.

ad_x = aprobe_df['rminrsep_D']
bd_x = bprobe_df['rminrsep_D']
cd_x = cprobe_df['rminrsep_D']
ad_y = aprobe_df['w_areal_D']
bd_y = bprobe_df['w_areal_D']
cd_y = cprobe_df['w_areal_D']
au_x = aprobe_df['rminrsep_U']
bu_x = bprobe_df['rminrsep_U']
cu_x = cprobe_df['rminrsep_U']
au_y = aprobe_df['w_areal_U']
bu_y = bprobe_df['w_areal_U']
cu_y = cprobe_df['w_areal_U']

def get_mean_std(xs, fs):
     min_idx = np.argmin(xs<10)
     max_idx = np.argmin(xs<15)
     mean = np.mean(fs[min_idx:max_idx])
     std  = np.std(fs[min_idx:max_idx])
     return mean, std

def get_other(frac):
     other_perp_lcfs = 0.5 * frac * cs_lcfs * ne_lcfs
     other_perp = other_perp_lcfs * np.exp(-r_range / lambda_gamma)
     return other_perp

ad_mean, ad_std = get_mean_std(x_ad, frac_ad)
bd_mean, bd_std = get_mean_std(x_bd, frac_bd)
cd_mean, cd_std = get_mean_std(x_cd, frac_cd)
au_mean, au_std = get_mean_std(x_au, frac_au)
bu_mean, bu_std = get_mean_std(x_bu, frac_bu)
cu_mean, cu_std = get_mean_std(x_cu, frac_cu)

ad_other = get_other(ad_mean)
bd_other = get_other(bd_mean)
cd_other = get_other(cd_mean)
au_other = get_other(au_mean)
bu_other = get_other(bu_mean)
cu_other = get_other(cu_mean)

def makeplot():
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    #ax1 = fig.add_subplot(111)
    #ax2 = fig.add_subplot(211)
    ax1.semilogy(ad_x, ad_y*1e19/125, 'r.', label='A-ITF')
    ax1.semilogy(bd_x, bd_y*1e19/125, 'b.', label='B-ITF')
    ax1.semilogy(cd_x, cd_y*1e19/125, 'g.', label='C-ITF')
    ax1.semilogy(r_range*100, ad_other, 'r--')
    ax1.semilogy(r_range*100, bd_other, 'b--')
    ax1.semilogy(r_range*100, cd_other, 'g--')
    #ax1.set_xlabel('R-Rsep (cm)')
    ax1.set_ylabel('Tungsten Flux (m-2 s-1)')
    ax1.legend()
    ax2.semilogy(au_x, au_y*1e19/125, 'r.', label='A-OTF')
    ax2.semilogy(bu_x, bu_y*1e19/125, 'b.', label='B-OTF')
    ax2.semilogy(cu_x, cu_y*1e19/125, 'g.', label='C-OTF')
    ax2.semilogy(r_range*100, au_other, 'r--')
    ax2.semilogy(r_range*100, bu_other, 'b--')
    ax2.semilogy(r_range*100, cu_other, 'g--')
    ax2.set_xlabel('R-Rsep (cm)')
    ax2.set_ylabel('Tungsten Flux (m-2 s-1)')
    ax2.legend()
    fig.tight_layout()
