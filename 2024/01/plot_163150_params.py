import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# Parameters from input file.
eV = 1.602e-19
qe, qi = -eV, eV
me, mp = 9.11e-31, 1.67e-27
AMU = 2.01410177811
mi = mp * AMU
mc = 12.011 * AMU
Te0 = 200 * eV
Ti0 = 200 * eV
n0 = 2.0e19
Zc = 6
fc = 0.02
qc = eV * Zc
Tc0 = Ti0
nc0 = n0 * fc
ni0 = n0 - nc0 * Zc
a_shift = 0.6
Z_axis = 0.00232616113
R_axis = 1.72068012
B_axis = 2.0
R_LCFSmid = 2.2801477223421736
Rmid_min = R_LCFSmid - 0.1
Rmid_max = R_LCFSmid + 0.05
R0 = 0.5 * (Rmid_min + Rmid_max)
a_mid = R_axis / a_shift - np.sqrt(R_axis * (R_axis - 2 * a_shift * R_LCFSmid + 2 * a_shift * R_axis)) / a_shift
r0 = R0 - R_axis
B0 = B_axis * (R_axis / R0)
kappa = 1.488
delta = 0.1


def r_x(x):
    """
    Go from simulation coordinate to minor radius.
    """
    return x + a_mid-0.1


def R(r, theta):
    return R_axis + r * np.cos(theta + np.arcsin(delta) * np.sin(theta))


def Z(r, theta):
    return Z_axis + kappa * r * np.sin(theta)


def Bphi(R):
    return B0 * R0 / R


def dZdr(r, theta):
    return kappa * np.sin(theta)


def dZdtheta(r, theta):
    return kappa * r * np.cos(theta)


def dRdr(r, theta):
    return np.cos(np.arcsin(delta) * np.sin(theta) + theta)


def dRdtheta(r, theta):
    return -r * np.sin(np.arcsin(delta) * np.sin(theta) + theta) * (np.arcsin(delta) * np.cos(theta) + 1)


def Jr(r, theta):
    # df is differentiation of the function. So df(R,1) is differentiation of R wrt to its first argument (r) and
    # df(R,2) is differentiation wrt to its second argument (theta). We just implement functions where I have calculated
    # the derivatives (e.g., WolframAlpha) and not rely on automatic differentiation.
    # return R(r,theta)*(df(R,1)(r,theta)*df(Z,2)(r,theta) - df(Z,1)(r,theta)*df(R,2)(r,theta))
    # print("dRdr = {}".format(dRdr(r, theta)))
    # print("dZdt = {}".format(dZdtheta(r, theta)))
    # print("dZdr = {}".format(dZdr(r, theta)))
    # print("dRdt = {}".format(dRdtheta(r, theta)))
    return R(r, theta) * (
            dRdr(r, theta) * dZdtheta(r, theta) - dZdr(r, theta) * dRdtheta(r, theta))


def dPsidr(r, theta):
    def integrand(t): return Jr(r, t) / R(r, t) ** 2

    # print(r)
    # print(theta)
    integral, _ = integrate.quad(integrand, 0, 2 * np.pi, epsrel=1e-10)
    # print(integral)
    # print(_)
    return B0 * R_axis / (2 * np.pi * qprofile(r)) * integral


# This one required some rewriting.
def alpha(r, theta, phi):
    def integrand(t):
        return Jr(r, t) / R(r, t) ** 2

    t = theta
    while t < -np.pi:
        t = t + 2 * np.pi
    while np.pi < t:
        t = t - 2 * np.pi
    if 0 < t:
        # integral, _ =  quad.dblexp(integrand, 0, t, 1e-10)/dPsidr(r,theta)
        integral, _ = integrate.quad(integrand, 0, t, epsabs=1e-10)
        integral = integral / dPsidr(r, theta)
    else:
        # integral, _ = -quad.dblexp(integrand, t, 0, 1e-10)/dPsidr(r,theta)
        integral, _ = integrate.quad(integrand, t, 0, epsabs=1e-10)
        integral = integral / dPsidr(r, theta)

    return phi - B0 * R_axis * integral


def gradr(r, theta):
    # return R(r,theta)/Jr(r,theta)*np.sqrt(df(R,2)(r,theta)^2 + df(Z,2)(r,theta)^2)
    return R(r, theta) / Jr(r, theta) * np.sqrt(
        dRdtheta(r, theta) ** 2 + dZdtheta(r, theta) ** 2)


def qprofile(r):
    a = [384.2108662176567, -2417.1050263480997, 5079.5143390131225, -3563.038325842346]
    return a[0] * (r + R_axis) ** 3 + a[1] * (r + R_axis) ** 2 + a[2] * (r + R_axis) + a[3]


def Bt(x, y, z):
    r = r_x(x)
    return Bphi(R(r, z))

def Bp(x, y, z):
    r = r_x(x)
    return dPsidr(r, z) / R(r, z) * gradr(r, z)

def bmag(x, y, z):
    r = r_x(x)
    Bt = Bphi(R(r, z))
    Bp = dPsidr(r, z) / R(r, z) * gradr(r, z)
    return np.sqrt(Bt ** 2 + Bp ** 2)


# More input file stuff
vti, vte, vtc = np.sqrt(Ti0 / mi), np.sqrt(Te0 / me), np.sqrt(Tc0 / mc)
c_s = np.sqrt(Te0 / mi)
omega_ci = np.abs(qi * B0 / mi)
rho_s = c_s / omega_ci
Lx = Rmid_max - Rmid_min
xMin, xMax = 0., Lx
rMin, rMax = Rmid_min - R_axis, Rmid_max - R_axis
q0 = qprofile(r_x(0.5 * (xMin + xMax)))
Ly = 150 * rho_s
Lz = 2. * np.pi + 1e-8
x_LCFS = R_LCFSmid - Rmid_min
ntoroidal = 2 * np.pi * r0 / q0 / Ly

# Create some profiles.
xs = np.linspace(xMin, xMax, 100)
zs = np.linspace(-Lz/2, Lz/2, 100)
q = qprofile(r_x(xs))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

ax1.axvline(x_LCFS, color="k", linestyle="--")
ax1.plot(xs, q)
ax1.set_xlabel("x (m)")
ax1.set_ylabel("q")

ax2.plot(xs, [dPsidr(r, None) for r in r_x(xs)])
ax2.set_xlabel("x (m)")
ax2.set_ylabel("dPsir")

ax3.plot(zs, Bp(0, None, zs))
ax3.set_xlabel("theta")
ax3.set_ylabel("bmag")

fig.tight_layout()
fig.show()