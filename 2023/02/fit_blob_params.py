# This is a bit crazy of thing for me, but the purpose is to extract the inner blob properties by fitting the measured
# radial velocity to Eq. 5.1.13 in Theiler's thesis.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import griddata, interp1d
import pickle
import openadas


# User-prescribed inputs.
shot = 167195
instab_ratio = 10
# yinst = 100
deltn_n = 1.0
temult = 1.0

# Load in files.
if shot == 167195:
    ca_path1 = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_167195_1.tab"
    ca_path2 = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_167195_2.tab"
    mafot_m1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp4500_conn_-1.dat"
    mafot_p1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp4500_conn_+1.dat"
    gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
elif shot == 190484:
    ca_path1 = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_190484_1.tab"
    ca_path2 = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_190484_2.tab"
    mafot_m1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190484/lam_rcp3000_conn_-1.dat"
    mafot_p1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190484/lam_rcp3000_conn_+1.dat"
    gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190484/190484_3000.pickle"

# Constants.
mi_ev = (2 * 931.49e6) / np.square(3e8)  # eV s2/m2
mi_kg = 2 * 1.66e-27
elec = 1.602e-19

# Load conditional averaging results.
ca1 = pd.read_csv(ca_path1, delimiter="\t")
ca2 = pd.read_csv(ca_path2, delimiter="\t")
ca_r = np.append(ca1["R(cm)"], ca2["R(cm)"]) / 100
ca_vr = np.append(ca1["Vr(m/s)"], ca2["Vr(m/s)"])
ca_a = np.append(ca1["D_rad(cm)"], ca2["D_rad(cm)"]) / 2.0 / 100  # diameter to radius, cm to m
ca_te = np.append(ca1["Te(eV)"], ca2["Te(eV)"]) * temult
ca_ne = np.append(ca1["Ne (e18m-3)"], ca2["Ne (e18m-3)"]) * 1e18

# Remove non-zero data...
mask = ca_vr > 0
ca_r = ca_r[mask]
ca_vr = ca_vr[mask]
ca_a = ca_a[mask]
ca_te = ca_te[mask]
ca_ne = ca_ne[mask]

# Load connection length data.
def load_mafot(itf_path, otf_path):
    print("Loading MAFOT data...")
    columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
               "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
    mafot_itf = pd.read_csv(itf_path, skiprows=52, names=columns, delimiter="\t")
    mafot_otf = pd.read_csv(otf_path, skiprows=52, names=columns, delimiter="\t")
    conns_r = mafot_itf["R (m)"]
    conns_l_itf = mafot_itf["Lconn (km)"].values * 1000  # km to m
    conns_l_otf = mafot_otf["Lconn (km)"].values * 1000
    return {"r":conns_r, "litf":conns_l_itf, "lotf":conns_l_otf}


# The equation uses the wall-to-wall value for the connection length. Get that and interpolate at the CA locations.
mafot = load_mafot(mafot_m1_path, mafot_p1_path)
conns = mafot["litf"] + mafot["lotf"]
f_conn = interp1d(mafot["r"], conns)
ca_l = f_conn(ca_r)

# Load gfile for B data.
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]
B = np.sqrt(np.square(gfile["Bt"]) + np.square(gfile["Bz"]))

# Get B field value at each RCP point.
ca_rzs = zip(ca_r, np.full(len(ca_r), -0.185))
ca_b = griddata((Rs.flatten(), Zs.flatten()), B.flatten(), list(ca_rzs))

class FitClass:

    def __init__(self, ab, B, L, te, vb, r, ne):
        self.ab = ab
        self.B = B
        self.L = L
        self.te = te
        self.vb = vb
        self.r = r
        self.ne = ne
    def gen_vr(self, r, nu_iz):
        """
        Equation 5.1.13 from Theiler thesis.
        """

        ab = self.ab
        B = self.B
        L = self.L
        te = self.te
        cs = np.sqrt(te / mi_ev)
        larmor = np.sqrt(te * mi_kg * elec) / (elec * B)
        self.larmor = larmor
        self.cs = cs
        numerator = np.sqrt(2 * ab / r) * cs * deltn_n
        denominator = instab_ratio + (1 / (np.square(larmor) * L)) * np.sqrt(r / 2) * np.power(ab, 5/2) \
                      + nu_iz * np.sqrt(r * ab) / (np.sqrt(2) * cs)
        return numerator / denominator

    def solve_nn(self):

        ab = self.ab
        r = self.r
        te = self.te
        B = self.B
        L = self.L
        vb = self.vb
        ne = self.ne
        cs = np.sqrt(te / mi_ev)
        larmor = np.sqrt(te * mi_kg * elec) / (elec * B)
        yint = np.sqrt(2 / (r * ab)) * cs
        yinst = vb / ab
        nu_iz = (np.sqrt(2 * ab / r) * cs * deltn_n / vb - (yinst / yint) - (1 / (np.square(larmor) * L))
                * np.sqrt(r / 2) * np.power(ab, 5/2)) * np.sqrt(2) * cs / np.sqrt(r * ab)

        oa = openadas.OpenADAS()
        rate_df = oa.read_rate_coef_unres("/Users/zamperini/My Drive/Research/Data/openadas/scd96_h.dat")
        ioniz_rate_coef = oa.get_rate_coef(rate_df, te, ne * deltn_n)
        nn = nu_iz / ioniz_rate_coef
        self.nu_iz = nu_iz
        self.ioniz_rate_coef = ioniz_rate_coef
        self.nn = nn
        self.yint = yint
        self.yinst = yinst
        return nn

# Perform curve fitting procedure.
fit = FitClass(ca_a, ca_b, ca_l, ca_te, ca_vr, ca_r, ca_ne)
# bounds=((1, 0), (100, 1000)),
popt, pcov = curve_fit(fit.gen_vr, ca_r, ca_vr, p0=(10), bounds=((0), (np.inf)), maxfev=100000)
fit_vr = fit.gen_vr(ca_r, *popt)

#print("Blob Te ~ {:.1f}".format(popt[0]))
#print("nu_iz ~ {:.2e}".format(popt[1]))

fig, ax = plt.subplots(figsize=(5, 4))

ax.plot([0, 1100], [0, 1100], color="k", linestyle="--", zorder=5)
ax.scatter(fit_vr, ca_vr, zorder=10)

fig.tight_layout()
fig.show()