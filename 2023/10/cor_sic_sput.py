import numpy as np
import matplotlib.pyplot as plt



# D --> Al2O3 yield (keV, Y)
d_al2o3_e = np.array([0.100,   0.170,   0.250,   0.300,   0.500,   1.000,   2.000,   4.000,   8.000])
d_al2o3_y = np.array([0.00240, 0.00500, 0.01950, 0.02300, 0.01900, 0.03890, 0.04700, 0.02860, 0.01400])

# D -- > SiC yield (keV, Y)
d_sic_e = np.array([0.1, 0.1, 0.12, 0.2, 0.2, 0.25, 0.3, 0.3, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3,
                    4, 5, 5, 8])
d_sic_y = np.array([0.0126, 0.013, 0.0155, 0.0198, 0.0198, 0.02280, 0.0320, 0.0288, 0.0299, 0.04, 0.0278, 0.0392,
                    0.0374, 0.04, 0.0392, 0.0496, 0.0452, 0.025, 0.042, 0.0284, 0.0475, 0.0172, 0.0211, 0.0171, 0.03,
                    0.0317, 0.0212, 0.0171, 0.0206, 0.0055])

def sput_yield(Eion, ion, target):
    """

    """

    if ion == "deuterium":
        m1 = 2.014
        z1 = 1
    else:
        m1 = 0
        z1 = 0

    if target.lower() in ["aluminum oxide ", "corundum", "al2o3"]:

        # As a compound the fits assumed their average atomic mass for M2 (and assuming Z2 as well).
        m2 = (26.982 * 2 + 15.999 * 3) / 5
        z2 = (13 * 2 + 8 * 3) / 5

        if ion == "deuterium":
            Eth = 66.0
            Q = 0.141
        else:
            Eth = 0.0
            Q = 0.0

    elif target.lower() in ["sic", "silicon carbide"]:
        m2 = (28.086 + 12.011) / 2
        z2 = (14 + 6) / 2

        if ion == "deuterium":
            Eth = 30.1
            Q = 0.119
        else:
            Eth = 0.0
            Q = 0.0

    else:
        m2 = 0
        z2 = 0
        Eth = 0.0
        Q = 0.0

    Etf = 30.74 * (m1 + m2) / m2 * z1 * z2 * (z1**(2/3) + z2**(2/3))**(1/2)
    eps = Eion / Etf
    sn = 0.5 * np.log(1 + 1.2288 * eps) / (eps + 0.1728 * np.sqrt(eps) + 0.008 * eps ** 0.1504)
    return Q * sn * (1 - (Eth / Eion)**(2/3)) * (1 - Eth / Eion)**2


def yield_al2o3(Eion, atom, ion="deuterium"):
    """
    Calculate yield of either Al or O from Al2O3 assuming stociometetry (for every 2 Al sputtered you sputter 3 O, i.e.,
    no surface enrichment occurs).
    """

    if atom == "al":
         # factor = 26.982 * 2 / (26.982 * 2 * 15.999 * 3)
         factor = 2 / (2 + 3)
    elif atom == "o":
        # factor = 15.999 * 2 / (26.982 * 2 * 15.999 * 3)
        factor = 3 / (2 + 3)
    else:
        return None

    return sput_yield(Eion, ion, "al2o3") * factor


def yield_sic(Eion, atom, ion="deuterium"):

    if atom == "si":
        # factor = 28.086 / (28.086 + 12.011)
        factor = 0.5
    elif atom == "c":
        # factor = 12.011 / (28.086 + 12.011)
        factor = 0.5
    else:
        return None
    return sput_yield(Eion, ion, "sic") * factor


# Validation plots.
if True:
    eions = np.geomspace(0.1, 20, 100)
    yield_al2o3 = [sput_yield(e * 1000, "deuterium", "al2o3") for e in eions]
    yield_sic = [sput_yield(e * 1000, "deuterium", "sic") for e in eions]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True, sharey=True)
    ax1.set_title("Al2O3")
    ax1.scatter(d_al2o3_e, d_al2o3_y, zorder=10, edgecolors="k")
    ax1.plot(eions, yield_al2o3, zorder=5)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid(which="both", zorder=1, alpha=0.3)
    ax1.set_xlabel("Energy (keV)")
    ax1.set_ylabel("Yield")
    ax2.set_title("SiC")
    ax2.scatter(d_sic_e, d_sic_y, zorder=10, edgecolors="k")
    ax2.plot(eions, yield_sic, zorder=5)
    ax2.grid(which="both", zorder=1, alpha=0.3)
    ax2.set_xlabel("Energy (keV")
    fig.tight_layout()
    fig.show()


