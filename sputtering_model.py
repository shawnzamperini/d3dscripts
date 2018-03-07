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
    #E0 = float(E0)

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
    elif ion == "carbon":
        # Carbon parameters are not readily avaible, so they were determined
        # using the function below (i.e. fitting to TRIM data).
        m1   = 12.01
        Z1   = 6
        q    = 5.4744
        mu   = 2.0321
        lamb = 20.5545
        Eth  = 30.0
    else:
        print("Error: Ion must be one of:\n  Deuterium\n  Tungsten\n  Carbon")
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
        for index in range(0, len(E0)):
            print("Energy = {} eV".format(E0[index]))
            print("  Lindhard screening length = {:.4e}".format(lindhard()[index]))
            print("  Reduced energy =            {:.4e}".format(reduced_energy()[index]))
            print("  Nuclear stopping power =    {:.4e}".format(nuclear_stopping()[index]))
            print("")

    # Return the yield using the above functions.
    ans = q * nuclear_stopping() * (E0 / Eth - 1)**mu / (lamb + (E0 / Eth - 1)**mu)
    if isinstance(ans, complex):
        return 0
    else:
        return ans


def Y_C_on_W(xmin=0, xmax=1000, plot_it1=False, plot_it2=False, return_it=True):
    # Constants needed.
    m1 = 12.01
    m2 = 183.84
    Z1 = 6
    Z2 = 74
    elec_sq = 1.44  # In eV*nm

    # Load csv file with pandas into a DataFrame, ignore the header info.
    yields = pd.read_csv('Data/C_on_W_yields.dat', skiprows=16, delim_whitespace=True)

    # Append column names to DataFrame.
    yields.columns = ['Energy', 'Alpha', 'Mean Depth', 'Refl. Coeff. C',
                      'Energy Coeff. C (2)', 'Sputt. Coeff. C', 'Energy Coeff. C (1)',
                      'Sputt. Coeff. W', 'Energy Coeff. W']

    # Run 1 is the first 18 indices. Run 2 is the rest.
    yields1 = yields[0:19]
    yields2 = yields[19:]

    def bohdansky(E0, q, Eth):
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

    def eckstein(E0, q, Eth, mu, lamb):
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
        #if isinstance(ans, complex):
        #    return 0
        #else:
        return ans

    # Our x-values are the energy, and y values are the yields.
    x1 = yields1['Energy'].values
    y1 = yields1['Sputt. Coeff. W'].values
    x1_fit = np.linspace(xmin, xmax, 1000)
    x2 = yields2['Energy'].values
    y2 = yields2['Sputt. Coeff. W'].values
    x2_fit = np.linspace(0, 400, 100)

    # The Bohdansky fit.
    popt1_boh, pcov1_boh = curve_fit(bohdansky, x1, y1)
    y1_fit_boh = bohdansky(x1_fit, *popt1_boh)
    popt2_boh, pcov2_boh = curve_fit(bohdansky, x2, y2)
    y2_fit_boh = bohdansky(x2_fit, *popt2_boh)

    # The Eckstein fit.
    guess = (5, 10, 2, 10)
    bounds = (-np.inf, [np.inf, np.inf, np.inf, np.inf])
    popt1_eck, pcov1_eck = curve_fit(eckstein, x1, y1, p0=guess, bounds=bounds)
    y1_fit_eck = eckstein(x1_fit, *popt1_eck)
    print("Fit #1")
    print("  q:    {:.2f}".format(popt1_eck[0]))
    print("  Eth:  {:.2f}".format(popt1_eck[1]))
    print("  mu:   {:.2f}".format(popt1_eck[2]))
    print("  lamb: {:.2f}".format(popt1_eck[3]))
    #guess = (5, 10, 2, 10)
    guess = popt1_eck
    bounds = (-np.inf, [np.inf, np.inf, np.inf, np.inf])
    popt2_eck, pcov2_eck = curve_fit(eckstein, x2[2:], y2[2:], p0=guess, maxfev=10000)
    y2_fit_eck = eckstein(x2_fit, *popt2_eck)
    print("Fit #2")
    print("  q:    {:.4f}".format(popt2_eck[0]))
    print("  Eth:  {:.4f}".format(popt2_eck[1]))
    print("  mu:   {:.4f}".format(popt2_eck[2]))
    print("  lamb: {:.4f}".format(popt2_eck[3]))

    if plot_it1:
        plt.plot(x1_fit, y1_fit_eck, 'r--', label="Eckstein")
        plt.plot(x1, y1, '.')
        plt.legend()
        plt.xlabel("Energy (eV)")
        plt.ylabel("Yield")
        plt.show()
    if plot_it2:
        plt.plot(x2_fit, y2_fit_eck, 'r--', label="Eckstein")
        #plt.plot(x2_fit, y2_fit_boh, 'b--', label="Bohdansky")
        plt.plot(x2, y2, '.')
        plt.legend()
        plt.xlabel("Energy (eV)")
        plt.ylabel("Yield")
        plt.show()
    if return_it:
        return x1_fit, y1_fit_eck, x2_fit, y2_fit_eck


def plot_of_yields(Emin=1, Emax=400):
    energies = np.linspace(Emin, Emax, 100)
    for ion in ["deuterium", "tungsten", "carbon"]:
        yields = Y(energies, ion=ion)
        plt.semilogy(energies, yields, label=ion)

    plt.xlabel("Energy (eV)")
    plt.ylabel("Yield")
    #plt.xlim([10, Emax])
    plt.legend()
    plt.show()

def yields(Emin=0.001, Emax=1000):
    energies = np.linspace(Emin, Emax, 10000)
    yields = {"deuterium": Y(energies, ion="deuterium"),
              "tungsten":  Y(energies, ion="tungsten"),
              "carbon":    Y(energies, ion="carbon")}
    y_frame = pd.DataFrame(yields, index=energies)

    # Anything that is blank should indicate zero yield.
    y_frame.fillna(0.0, inplace=True)
    y_frame.index.name = "Energy (eV)"
    return y_frame

def get_yield(y_frame, ion, charge_state=1, kTe=25, verbose=True):
    """
    Will return the yield of an ion at a specific charge state assuming Ti=Te.
    y_frame: DataFrame returned from "yields" function.
    """
    # Constants
    boltz = 8.617e-5 # Boltzmann constant in eV/K

    if ion == "carbon":
        if charge_state > 6:
            print("Error: Carbon charge state cannot be greater than 6.")
            return None
    if ion == "tungsten":
        if charge_state > 6:
            print("Error: Not considering tungsten above W6+.")

    E_impact = 2 * kTe + 3 * charge_state * kTe

    # Find closest energy in frame to E_impact.
    ion_series     = y_frame[ion]
    energies       = ion_series.index.values
    idx            = np.abs(energies - E_impact).argmin()
    closest_energy = energies[idx]
    if verbose:
        print("Impact energy is:  {0:.4f} eV".format(E_impact))
        print("Closest energy is: {0:.4f} eV".format(closest_energy))
        print("  Yield is {0:.4f}".format(ion_series[closest_energy]))

    return ion_series[closest_energy]
