import numpy as np


e = 1.609 * 10**(-19)
mw = 183.84 * 1.66 * 10**(-27)
q = 5 * e
Bt = 1.5
B = 1.5
E0 = 1500
dt = 2 * 10**(-5)
vr0 = 100

#dvr = (1/2) * (q / mw)**2 * (Bt*E0 - B**2*vr0) * dt**2

#dvr = Bt * E0 * np.sqrt(2 * np.pi) / np.square(B)
dvr = Bt * E0 / np.square(B)
dr = dvr * dt

print("dr = {:.3e}".format(dr))
