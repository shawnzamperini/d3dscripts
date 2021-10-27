# Simple thing to go from UOB volts to particles per second.
import numpy as np


uob_volts = 1.2  # Same for 184267 and 184527.

# Calibration constants.
p1 = 10.31
p2 = 0.38
flow3v = 5.26

# Flow in Torr-L/s.
f_tlps = p1 * (np.sqrt(np.square(uob_volts * p2) + 1) - 1)

# Flow in sccm.
f_sccm = f_tlps * 78.9

# Flow in particles/s.
f_pps = f_sccm * 4.477962 * 10**17
print("{:.1f} V".format(uob_volts))
print("  Way #1: {:.3e} atoms/s".format(f_pps))

# The ideal gas law way, PV = nRT
R = 62.363  # Ideal gas constant in L*Torr/K*mol.
T = 293  # Room temperature in Kelvin (293) or wall temperature (373).
n = f_tlps / (R * T)  # Number of moles/s.
avog = 6.02214 * 10**(23)  # Avogadro's constant.
N = n * avog  # Number of molecules/s.
print("  Way #2: {:.3e} atoms/s".format(N))
