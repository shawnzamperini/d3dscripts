import numpy as np

amu = 1.660e-27  # kg
mi_amu = 2.014
mi = amu * mi_amu  # kg
#mW_amu = 183.38
mW_amu = mi_amu
#mW = amu * mW_amu
mW = mi
me = 9.10938188e-31  # kg
ev = 1.602e-19  # C
lnalpha = 15
eps0 = 8.85e-12

#ZW = 5 # W charge state
ZW = 1
Z = 1  # Deuterium
Ti = 10
Te = 10
ni = 1e18
ne = ni
B = 1.5  # T

vtW = np.sqrt(2 * Ti * ev / mW)
vti = np.sqrt(2 * Ti * ev / mi)
vte = np.sqrt(2 * Te * ev / me)
mr = mi * mW / (mi + mW)
mr_amu = mi_amu * mW_amu / (mi_amu + mW_amu)
mr_ei = mi * me / (mi + me)
mr_ii_amu = mi_amu * mi_amu / (mi_amu + mi_amu)

# From Helander & Sigmar
tau_zi = (12 * np.power(np.pi, 3/2) * np.sqrt(mW) * np.power(Ti*ev, 3/2) * np.square(eps0)) / (np.sqrt(2) * ni * np.square(ZW) * np.square(Z) * np.power(ev, 4) * lnalpha)

# From Friedberg, Eq. 9.48. Only valid for W1+, missing correction.
nu_fr = np.power(ev, 4) * ni * lnalpha / (4 * np.pi * np.square(eps0) * mW * mr) / (np.power(vtW, 3) + 1.3 * np.power(vti, 3))

# From Friedberg again, but with impurity charge state and corrected lnalpha.
lnalpha_corr = np.log(np.sqrt(eps0 * Te * ev / (ne * np.square(ev))) * 4 * np.pi 
    * eps0 * mr * np.square(vtW) / (ZW * Z * np.square(ev)))
lnalpha_fact = np.log(np.power(eps0, 3/2) * 4 * np.pi * amu / np.power(ev, 5/2))
print("lnalpha_fact = {:}".format(lnalpha_fact))
lnalpha_corr2 = lnalpha_fact + np.log(np.sqrt(Te / ne) * np.square(vtW) * mr_amu / (ZW * Z))
print("lnalpha_corr = {:.2f}".format(lnalpha_corr))
print("lnalpha_corr2 = {:.2f}".format(lnalpha_corr2))
lnalpha_ei = np.log(np.sqrt(eps0 * Te * ev / (ne * np.square(ev))) * 4 * np.pi 
    * eps0 * mr_ei * np.square(vte) / (np.square(ev)))
print("lnalpha_ei = {:.2f}".format(lnalpha_ei))
nu_fr_corr = np.square(ZW*ev) * np.square(Z*ev) * ni * lnalpha_corr / (4 * np.pi 
    * np.square(eps0) * mW * mr) / (np.power(vtW, 3) + 1.3 * np.power(vti, 3))
nu_fr_fact = np.power(ev, 4) * np.power(amu, 3/2) / (4 * np.pi * np.square(eps0) * np.square(amu) * np.power(2 * ev, 3/2))
print("nu_fr_fact = {:.2e}".format(nu_fr_fact))
nu_fr_corr2 = nu_fr_fact * np.square(ZW) * np.square(Z) * ni * lnalpha_corr2 / (mW_amu * mr_amu) / (np.power(Ti / mW_amu, 3/2) + 1.3 * np.power(Ti / mi_amu, 3/2))
print("nu_fr_corr2 = {:.2e}".format(nu_fr_corr2))

# The ion-ion collision frequency for comparison.
lnalpha_corr_ii = lnalpha_fact + np.log(np.sqrt(Te / ne) * np.square(vti) * mr_amu / (Z * Z))
nu_fr_corr_ii = nu_fr_fact * np.square(Z) * np.square(Z) * ni * lnalpha_corr_ii / (mi_amu * mr_ii_amu) / (np.power(Ti / mi_amu, 3/2) + 1.3 * np.power(Ti / mi_amu, 3/2))

print("H&S")
print("  tau_zi = {:.3e}".format(tau_zi))
print("  nu_zi = {:.3e}".format(1 / tau_zi))
print("Friedberg")
print("  tau_zi = {:.3e}".format(1 / nu_fr))
print("  nu_zi = {:.3e}".format(nu_fr))
print("Friedberg - Corrected")
print("  tau_zi = {:.3e}".format(1 / nu_fr_corr))
print("  nu_zi = {:.3e}".format(nu_fr_corr))
print("  nu_ii = {:.3e}".format(nu_fr_corr_ii))
print("lambda_mfpW = {:.3e}".format(1 / nu_fr_corr2 * vtW))
rhoW = vtW / (ZW * ev * B / mW)
print("vtW = {:.2e}".format(vtW))
print("rhoW = {:.2e}".format(rhoW))

