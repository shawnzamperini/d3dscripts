import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def Y(E0, ion, debug=False):
    """
    Fitted yield from Eckstein, assuming normal incidence. User only needs to
    enter the incident ion energy, E0, and what the ion is, either "deuterium"
    or "carbon".
        E0:
            Incident ion energy (eV).
        m1, m2:
            Mass of projectile and target atom, respectively.
        Z1, Z2:
            Atomic number of projectile and target atom, respectively.
        Eth, q, mu, lamb:
            Fitting parameters determined from Eckstein.
            Eth: Threshold energy for sputtering (eV).
            q: Related to absolute yield.
            lamb: Related to onset of decrease of yield at low energies towards
                  the threshold.
            mu: Describes strength of this decrease.
    """

    # Make sure in lowercase.
    ion = ion.lower()
    # Some weird error only with numpy.float64. Just use floats I guess.
    E0 = float(E0)

    # Masses in amu since all it uses it the ratio between them.
    m2 = 183.84
    Z2 = 74
    elec_sq = 1.44  # In eV*nm

    # Select the relevant fitting parameters.
    if ion == "deuterium":
        m1   = 2.014
        Z1   = 1
        q    = 0.0183
        mu   = 1.4410
        lamb = 0.3583
        Eth  = 228.84
    elif ion == "tungsten":
        m1   = 183.84
        Z1   = 74
        q    = 18.6006
        mu   = 3.1273
        lamb = 2.2697
        Eth  = 24.9885
    else:
        print("Error: Ion must be one of:\n  Deuterium\n  Tungsten")
        return None


    def lindhard():
        """ The Lindhard screening length (nm)."""
        return (9 * np.pi**2 / 128)**(1/3) * 0.0529177 * (Z1**(2/3) + Z2**(2/3))**(-1/2)

    def reduced_energy():
        """ The reduced energy (unitless). """
        return E0 * m1/(m1+m2) * lindhard() / (Z1 * Z2 * elec_sq)

    def nuclear_stopping():
        """ The nuclear stopping power with KrC (WHB) potential (). """
        return 0.5 * np.log(1 + 1.2288 * reduced_energy()) / (reduced_energy() +
               0.1728 * np.sqrt(reduced_energy()) + 0.008 * reduced_energy()**0.1504)

    if debug:
        print("Energy = {} eV".format(E0))
        print("  Lindhard screening length = {:.4e}".format(lindhard()))
        print("  Reduced energy =            {:.4e}".format(reduced_energy()))
        print("  Nuclear stopping power =    {:.4e}".format(nuclear_stopping()))
        print("")

    # Return the yield using the above functions.
    ans = q * nuclear_stopping() * (E0 / Eth - 1)**mu / (lamb + (E0 / Eth - 1)**mu)
    if isinstance(ans, complex):
        return 0.0
    else:
        return ans


def Y_C_on_W():

    # Load csv file with pandas into a DataFrame, ignore the header info.
    yields = pd.read_csv('Data/C_on_W_yields.dat', skiprows=16, delim_whitespace=True)

    # Append column names to DataFrame.
    yields.columns = ['Energy', 'Alpha', 'Mean Depth', 'Refl. Coeff. C',
                      'Energy Coeff. C (2)', 'Sputt. Coeff. C', 'Energy Coeff. C (1)',
                      'Sputt. Coeff. W', 'Energy Coeff. W']

    # Run 1 is the first 18 indices. Run 2 is the rest.
    yields1 = yields[0:19]
    yields2 = yields[19:]

    # Want to fit the E vs. Y with Eckstein fitting and compare.
    def eckstein(E0, q, mu, lamb, Eth):
        # Constants needed.
        m1 = 12.01
        m2 = 183.84
        Z1 = 6
        Z2 = 74
        elec_sq = 1.44  # In eV*nm

        def lindhard():
            """ The Lindhard screening length (nm)."""
            return (9 * np.pi**2 / 128)**(1/3) * 0.0529177 * (Z1**(2/3) + Z2**(2/3))**(-1/2)
        def reduced_energy():
            """ The reduced energy (unitless). """
            return E0 * m1/(m1+m2) * lindhard() / (Z1 * Z2 * elec_sq)
        def nuclear_stopping():
            """ The nuclear stopping power with KrC (WHB) potential (). """
            return 0.5 * np.log(1 + 1.2288 * reduced_energy()) / (reduced_energy() +
                   0.1728 * np.sqrt(reduced_energy()) + 0.008 * reduced_energy()**0.1504)
        # Return the yield using the above functions.
        ans = q * nuclear_stopping() * (E0 / Eth - 1)**mu / (lamb + (E0 / Eth - 1)**mu)

        if False:
            print("Energy = {} eV".format(E0))
            print("  Lindhard screening length = {}".format(lindhard()))
            print("  Reduced energy =            {}".format(reduced_energy()))
            print("  Nuclear stopping power =    {}".format(nuclear_stopping()))
            print("")

        if isinstance(ans, complex):
            return 0.0
        else:
            return ans

    def bohdansky(E0, q, Eth):
        # Constants needed.
        m1 = 12.01
        m2 = 183.84
        Z1 = 6
        Z2 = 74
        elec_sq = 1.44  # In eV*nm
        def lindhard():
            """ The Lindhard screening length (nm)."""
            return (9 * np.pi**2 / 128)**(1/3) * 0.0529177 * (Z1**(2/3) + Z2**(2/3))**(-1/2)
        def reduced_energy():
            """ The reduced energy (unitless). """
            return E0 * m1/(m1+m2) * lindhard() / (Z1 * Z2 * elec_sq)
        def nuclear_stopping():
            """ The nuclear stopping power with KrC (WHB) potential (). """
            return 0.5 * np.log(1 + 1.2288 * reduced_energy()) / (reduced_energy() +
                   0.1728 * np.sqrt(reduced_energy()) + 0.008 * reduced_energy()**0.1504)
        ans = q * nuclear_stopping() * (1 - Eth/E0**(2/3)) * (1 - Eth/E0)**2
        if isinstance(ans, complex):
            return 0.0
        else:
            return ans

    # Our x-values are the energy, and y values are the yields.
    x1 = yields1['Energy'].values
    y1 = yields1['Sputt. Coeff. W'].values
    if True:
        guess = (10, 2, 0.1, 50)
        popt, pcov = curve_fit(eckstein, x1, y1, p0=guess, maxfev=50000, bounds=(0, [np.inf, np.inf, 10, np.inf]))
        print("q, mu, lambda, Eth: {0}, {1}, {2}, {3}".format(popt[0], popt[1], popt[2], popt[3]))
        y1_fit = eckstein(x1, *popt)
        plt.plot(x1, y1_fit, 'b--', label="Eckstein")
    if True:
        popt, pcov = curve_fit(bohdansky, x1, y1)
        print("q, Eth: {0}, {1}".format(popt[0], popt[1]))
        y1_fit = bohdansky(x1, *popt)
        plt.plot(x1, y1_fit, 'r--', label="Bohdansky")
    plt.plot(x1, y1, '.')
    plt.legend()
    plt.xlabel("Energy (eV)")
    plt.ylabel("Yield")
    plt.show()


def plot_of_yields(ion, Emin=1, Emax=10000):
        energies = np.linspace(Emin, Emax, 1000)
        yields   = np.array([Y(energy, ion) for energy in energies])

        plt.loglog(energies, yields)

        if ion.lower() == "deuterium":
            title = "D on W Yield"
        elif ion.lower() == "tungsten":
            title = "W on W Yield"
        plt.title(title)

        plt.xlabel("Energy (eV)")
        plt.ylabel("Yield")

        # Only plot over the range when Y > 0.
        low_index = np.argmax(yields > 0)
        plt.xlim(energies[low_index], np.max(energies))

        plt.show()
