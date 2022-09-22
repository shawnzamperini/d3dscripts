import numpy as np


# Inputs.
dperp = 1.0
L = 10.0
te = 5
ti = te
x0 = 1
rdpout = 28  # Torr-L/s, jumps from ~28 to ~32

# Derived quantities.

# The ideal gas law way, PV = nRT
R = 62.363  # Ideal gas constant in L*Torr/K*mol.
T = 293  # Room temperature in Kelvin (293) or wall temperature (373).
n = rdpout / (R * T)  # Number of moles/s.
avog = 6.02214 * 10**(23)  # Avogadro's constant.
S0 = n * avog  # Number of molecules/s.
print("  S0 = {:.3e} atoms/s".format(S0))

mi = 931.49e6
cs = np.sqrt((te+ti)/mi) * 3e8


print("Term #1: {:.2e}".format(np.sqrt(dperp*L/cs)))
print("Term #2: {:.2e}".format(1/np.sqrt((np.pi-2)/2 + S0*x0/cs)))
lambdane = np.sqrt(dperp*L/cs) * 1 / np.sqrt((np.pi-2)/2 + S0*x0/cs)
