import sys
sys.path.append("/home/zamperinis/d3dscripts/Utilities")
import openadas
import numpy as np


# Create OpenADAS object 
oa = openadas.OpenADAS()

# DataFrames for ionization (scd) and recombination (acd) rates
w_scd = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/scd50/scd50_w.dat")
w_acd = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/acd50/acd50_w.dat")

# Radiation values calculated for a typical core plasma
core_te = 5e3
core_ne = 1e20 # Is this typical?

# Calculate fractional abundances for each impurity at core values
w_frac_abund = oa.coronal_equil_frac_abundance(74, w_scd, w_acd, core_te, core_ne)

# DataFrames for line radiation (plt) and bremstrahlung (prb)
w_plt = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/plt50/plt50_w.dat")
w_prb = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/prb50/prb50_w.dat")
fe_plt = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/plt89/plt89_fe.dat")
fe_prb = oa.read_rate_coef_unres("/fusion/projects/adas/adas/adf11/prb89/prb89_fe.dat")

# Get the respective powers. I can only assume these are W m3 (they've been internally
# multiplied within the OpenADAS class by 1e-6 to do cm3 --> m3). 
w_plt_rates = np.array([oa.get_rate_coef(w_plt, core_te, core_ne, charge=z+1) for z in range(74)])
w_prb_rates = np.array([oa.get_rate_coef(w_prb, core_te, core_ne, charge=z+1) for z in range(74)])

# Example concentration
w_conc = 1e-3

# To calculate W/m3 it's rate[Z] * ne * nz[Z] = rate[Z] * ne * (ne * w_conc * frac[Z])
w_plt_rad = w_plt_rates * core_ne ** 2 * w_conc * w_frac_abund[1:]  # [:1] because skip neutrals
w_prb_rad = w_prb_rates * core_ne ** 2 * w_conc * w_frac_abund[1:]  # [:1] because skip neutrals

# Total is just sum of these two contributions
w_tot_rad = w_plt_rad + w_prb_rad
