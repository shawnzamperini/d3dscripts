import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal   import savgol_filter
import warnings


def fit_to_eck(plot_it=True):
    # Constants needed.
    m1 = 12.01
    m2 = 183.84
    Z1 = 6
    Z2 = 74
    elec_sq = 1.44  # In eV*nm

    # Load csv file with pandas into a DataFrame, ignore the header info.
    yields_90C10W = pd.read_csv('../Data/C-CW_90C10W_sdtrimsp.dat', skiprows=16, delim_whitespace=True)
    yields_99C1W  = pd.read_csv('../Data/C-CW_99C1W_sdtrimsp.dat',  skiprows=16, delim_whitespace=True)
    yields_D_on_90C10W = pd.read_csv('../Data/D-CW_90C10W_sdtrimsp.dat', skiprows=16, delim_whitespace=True)
    yields_W_on_90C10W = pd.read_csv('../Data/W-CW_90C10W_sdtrimsp.dat', skiprows=16, delim_whitespace=True)

    # Append column names to DataFrame.
    col_names = ['Energy', 'Alpha', 'Mean Depth', 'Refl. Coeff. C', 'Energy Coeff. C (2)',
                 'Sputt. Coeff. C', 'Energy Coeff. C (1)', 'Sputt. Coeff. W', 'Energy Coeff. W']
    yields_90C10W.columns = col_names
    yields_99C1W.columns = col_names
    yields_D_on_90C10W.columns = ['Energy', 'Alpha', 'Mean Depth', 'Refl. Coeff. D', 'Energy Coeff. D (2)',
                                  'Sputt. Coeff. D', 'Energy Coeff. D (1)', 'Sputt. Coeff. C', 'Energy Coeff. C (1)',
                                  'Sputt. Coeff. W', 'Energy Coeff. W']
    yields_W_on_90C10W.columns = ['Energy', 'Alpha', 'Mean Depth', 'Refl. Coeff. W', 'Energy Coeff. W (2)',
                                  'Sputt. Coeff. W', 'Energy Coeff. W (1)', 'Sputt. Coeff. C', 'Energy Coeff. C']

    # The Eckstein equation to fit data to.
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

    all_yields = [yields_90C10W, yields_99C1W, yields_D_on_90C10W, yields_W_on_90C10W]
    for yields in all_yields:
        # Our x-values are the energy, and y values are the yields.
        x     = yields['Energy'].values
        y     = yields['Sputt. Coeff. W'].values
        x_fit = np.linspace(0, 500, 500)

        # The Eckstein fit.
        guess = (5, 10, 2, 10)
        bounds = (-np.inf, [np.inf, np.inf, np.inf, np.inf])
        popt, pcov = curve_fit(eckstein, x, y, p0=guess, bounds=bounds)
        y_fit_eck = eckstein(x_fit, *popt)
        if yields is yields_90C10W:
            print("Fit (90% C 10% W)")
            plt.title("90% C 10% W")
        elif yields is yields_99C1W:
            print("Fit (99% C 1% W)")
            plt.title("99% C 1% W")
        elif yields is yields_D_on_90C10W:
            print("Fit (D on 90% C 10%W)")
            plt.title("D on 90% C 10%W")
        elif yields is yields_W_on_90C10W:
            print("Fit (W on 90% C 10%W)")
            plt.title("W on 90% C 10%W")
        print("  q:    {:.2f}".format(popt[0]))
        print("  Eth:  {:.2f}".format(popt[1]))
        print("  mu:   {:.2f}".format(popt[2]))
        print("  lamb: {:.2f}".format(popt[3]))

        if plot_it:
            plt.plot(x_fit, y_fit_eck, 'r--', label="Eckstein")
            plt.plot(x, y, '.')
            plt.legend()
            plt.xlabel("Energy (eV)")
            plt.ylabel("Yield")
            plt.show()

    return all_yields

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
    #ion = ion.lower()
    # Some weird error only with numpy.float64. Just use floats I guess.
    #E0 = [float(E) for E in E0]

    # Masses in amu since all it uses it the ratio between them.
    m2 = 183.84
    Z2 = 74
    elec_sq = 1.44  # In eV*nm

    # Select the relevant fitting parameters.
    if ion == "deuterium_on_100W":
        m1   = 2.014
        Z1   = 1
        q    = 0.0183
        mu   = 1.4410
        lamb = 0.3583
        Eth  = 228.84
    elif ion == "tungsten_on_100W":
        m1   = 183.84
        Z1   = 74
        q    = 18.6006
        mu   = 3.1273
        lamb = 2.2697
        Eth  = 24.9885
    elif ion == "carbon_on_100W":
        # Carbon parameters are not readily available, so they were determined
        # using the function below (i.e. fitting to TRIM data).
        m1   = 12.01
        Z1   = 6
        q    = 5.4744
        mu   = 2.0321
        lamb = 20.5545
        Eth  = 30.0
    elif ion == "carbon_on_90C10W":
        m1   = 12.01
        Z1   = 6
        q    = 1.73
        mu   = 0.77
        lamb = 118.61
        Eth  = 50.0
    elif ion == "deuterium_on_90C10W":
        m1   = 2.014
        Z1   = 1
        q    = 0.01
        mu   = 1.90
        lamb = 87.41
        Eth  = 50.0
    elif ion == "tungsten_on_90C10W":
        m1   = 183.84
        Z1   = 74
        q    = -0.25
        mu   = 1.47
        lamb = -107.87
        Eth  = 75.0
    elif ion == "carbon_on_99C1W":
        m1   = 12.01
        Z1   = 6
        q    = 0.01
        mu   = 1.39
        lamb = 85.18
        Eth  = 50.0

    else:
        print("Error: Incorrect ion entry.")
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
    # Weird warning only involving np.floats, but doesn't seem to cause problems so ignore it.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ans = q * nuclear_stopping() * (E0 / Eth - 1)**mu / (lamb + (E0 / Eth - 1)**mu)
        if isinstance(ans, complex):
            return 0
        else:
            return ans

def yields(Emin=0.001, Emax=500):
    energies = np.linspace(Emin, Emax, 10000)
    yields = {"carbon_on_100W":        Y(energies, ion="carbon_on_100W"),
              "carbon_on_90C10W":      Y(energies, ion="carbon_on_90C10W"),
              "carbon_on_99C1W":       Y(energies, ion="carbon_on_99C1W"),
              "deuterium_on_90C10W":   Y(energies, ion="deuterium_on_90C10W"),
              "tungsten_on_90C10W":    Y(energies, ion="tungsten_on_90C10W")}
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

def get_fluxes(carbon_frac=  np.array([0.00,  0.005,  0.005,  0.00,  0.00,  0.00]),
               tungsten_frac=np.array([0.00005, 0.00005, 0.00, 0.00, 0.00, 0.00])):
               #carbon_frac):
               #tungsten_frac):
    """
    Calculates the fluxes of each ion charge state. The elements in the lists
    passed in each correspond to a charge state (i.e. C1+, C2+, ..., C6+).
    """

    excel_file = "/home/shawn/d3dscripts/Data/LP_with_fit.xlsx"

    if False:
        # Or use plunge 1 of 167192 since the separatrix location is more representative
        # of the separatrix during the cp shots (196-220).
        lp_df = pd.read_excel(excel_file,
                              sheet_name="LP Data",
                              skiprows=[0, 1],
                              names=["Time (ms)", "R (cm)", "ne (e18 m-3)", "Te (eV)",
                                     "Vfl (V)", "Vp (V)", "R-Rsep (cm)"],
                              usecols=[0, 1, 2, 3, 4, 5, 6])

        # Drop the null data it pulls.
        lp_df.drop(lp_df.index[56:], inplace=True)

    # Try using all the data (excludes sweeping plunges).
    if True:
        lp_df = pd.read_excel(excel_file,
                              sheet_name="Fit Te's",
                              #skiprows=[0],
                              names=['R (cm)', 'Te (eV)', 'ne (e18 m-3)'],
                              usecols=[0, 1, 3])

    # Want the flow at probe face (= flow at sheath edge). It's 0.5ne*cs
    m_deut      = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
    cs_series   = (2 * lp_df["Te (eV)"] / m_deut) ** (1/2)
    flux_series = 0.5 * lp_df['ne (e18 m-3)'] * cs_series * 10**(18)
    lp_df['D Flux (m-2 s-1)'] = flux_series

    # Estimate of the carbon and tungsten flux as fractions of the deuterium flux.
    for charge_state in range(0, 6):
        lp_df['C' + str(charge_state+1) +  '+ Flux (m-2 s-1)'] = \
            lp_df['D Flux (m-2 s-1)'] * carbon_frac[charge_state]
        lp_df['W' + str(charge_state+1) + '+ Flux (m-2 s-1)'] = \
            lp_df['D Flux (m-2 s-1)'] * tungsten_frac[charge_state]

    lp_df['R-Rsep (cm)'] = lp_df['R (cm)'] - 221.83

    return lp_df

def get_sput_flux(y_df, fluxes_df):
    #y_df = yields()
    #fluxes_df = get_fluxes()
    sput_df = pd.DataFrame()
    sput_df["R (cm)"] = fluxes_df["R (cm)"]
    sput_df["R-Rsep (cm)"] = sput_df['R (cm)'] - 221.83

    # Calculate the sputtered flux from each charge state of the ions.
    for ion in ["carbon_on_100W", "carbon_on_90C10W", "carbon_on_99C1W", "deuterium_on_90C10W", "tungsten_on_90C10W"]:
        # As of now only considering charge states 1-6.
        for charge_state in range(1, 7):
            if ion == "carbon_on_100W":
                fluxes = fluxes_df["C" +str(charge_state) + "+ Flux (m-2 s-1)"].values
                Tes    = fluxes_df["Te (eV)"].values
                sput_flux = np.array([])
                for flux, Te in zip(fluxes, Tes):
                    tmp_y    = get_yield(y_df, ion, charge_state, Te, verbose=False)
                    tmp_sput = tmp_y * flux
                    sput_flux = np.append(sput_flux, tmp_sput)
                sput_df['Sputt. Flux from C on 100W' + str(charge_state) + '+'] = sput_flux
            elif ion == "carbon_on_90C10W":
                fluxes = fluxes_df["C" +str(charge_state) + "+ Flux (m-2 s-1)"].values
                Tes    = fluxes_df["Te (eV)"].values
                sput_flux = np.array([])
                for flux, Te in zip(fluxes, Tes):
                    tmp_y    = get_yield(y_df, ion, charge_state, Te, verbose=False)
                    tmp_sput = tmp_y * flux
                    sput_flux = np.append(sput_flux, tmp_sput)
                sput_df['Sputt. Flux from C on 90C10W' + str(charge_state) + '+'] = sput_flux
            elif ion == "carbon_on_99C1W":
                fluxes = fluxes_df["C" +str(charge_state) + "+ Flux (m-2 s-1)"].values
                Tes    = fluxes_df["Te (eV)"].values
                sput_flux = np.array([])
                for flux, Te in zip(fluxes, Tes):
                    tmp_y    = get_yield(y_df, ion, charge_state, Te, verbose=False)
                    tmp_sput = tmp_y * flux
                    sput_flux = np.append(sput_flux, tmp_sput)
                sput_df['Sputt. Flux from C on 99C1W' + str(charge_state) + '+'] = sput_flux
            elif ion == "deuterium_on_90C10W":
                fluxes = fluxes_df["D Flux (m-2 s-1)"].values
                Tes    = fluxes_df["Te (eV)"].values
                sput_flux = np.array([])
                for flux, Te in zip(fluxes, Tes):
                    tmp_y    = get_yield(y_df, ion, 1, Te, verbose=False)
                    tmp_sput = tmp_y * flux
                    sput_flux = np.append(sput_flux, tmp_sput)
                sput_df['Sputt. Flux from D on 90C10W'] = sput_flux
            elif ion == "tungsten_on_90C10W":
                fluxes = fluxes_df["W" +str(charge_state) + "+ Flux (m-2 s-1)"].values
                Tes    = fluxes_df["Te (eV)"].values
                sput_flux = np.array([])
                for flux, Te in zip(fluxes, Tes):
                    tmp_y    = get_yield(y_df, ion, charge_state, Te, verbose=False)
                    tmp_sput = tmp_y * flux
                    sput_flux = np.append(sput_flux, tmp_sput)
                sput_df['Sputt. Flux from W on 90C10W' + str(charge_state) + '+'] = sput_flux

    sput_df["Sputt. Flux from C on 100W"] = sput_df["Sputt. Flux from C on 100W1+"].values + \
                                            sput_df["Sputt. Flux from C on 100W2+"].values + \
                                            sput_df["Sputt. Flux from C on 100W3+"].values + \
                                            sput_df["Sputt. Flux from C on 100W4+"].values + \
                                            sput_df["Sputt. Flux from C on 100W5+"].values + \
                                            sput_df["Sputt. Flux from C on 100W6+"].values
    sput_df["Sputt. Flux from C on 90C10W"] = sput_df["Sputt. Flux from C on 90C10W1+"].values + \
                                              sput_df["Sputt. Flux from C on 90C10W2+"].values + \
                                              sput_df["Sputt. Flux from C on 90C10W3+"].values + \
                                              sput_df["Sputt. Flux from C on 90C10W4+"].values + \
                                              sput_df["Sputt. Flux from C on 90C10W5+"].values + \
                                              sput_df["Sputt. Flux from C on 90C10W6+"].values
    sput_df["Sputt. Flux from C on 99C1W"]  = sput_df["Sputt. Flux from C on 99C1W1+"].values + \
                                              sput_df["Sputt. Flux from C on 99C1W2+"].values + \
                                              sput_df["Sputt. Flux from C on 99C1W3+"].values + \
                                              sput_df["Sputt. Flux from C on 99C1W4+"].values + \
                                              sput_df["Sputt. Flux from C on 99C1W5+"].values + \
                                              sput_df["Sputt. Flux from C on 99C1W6+"].values
    sput_df["Sputt. Flux from W on 90C10W"] = sput_df["Sputt. Flux from W on 90C10W1+"].values + \
                                              sput_df["Sputt. Flux from W on 90C10W2+"].values + \
                                              sput_df["Sputt. Flux from W on 90C10W3+"].values + \
                                              sput_df["Sputt. Flux from W on 90C10W4+"].values + \
                                              sput_df["Sputt. Flux from W on 90C10W5+"].values + \
                                              sput_df["Sputt. Flux from W on 90C10W6+"].values
    # Convert to fluences.
    sput_df["Sputt. Fluence from C on 100W"]   = sput_df["Sputt. Flux from C on 100W"]   * 125.0
    sput_df["Sputt. Fluence from C on 90C10W"] = sput_df["Sputt. Flux from C on 90C10W"] * 125.0
    sput_df["Sputt. Fluence from C on 99C1W"]  = sput_df["Sputt. Flux from C on 99C1W"]  * 125.0
    sput_df["Sputt. Fluence from W on 90C10W"] = sput_df["Sputt. Flux from W on 90C10W"] * 125.0
    sput_df["Sputt. Fluence from D on 90C10W"] = sput_df["Sputt. Flux from D on 90C10W"] * 125.0

    return sput_df

def calc_net_fluence():
    exposed_time = 5.0 * 25.0
    filename = '/home/shawn/d3dscripts/Data/LModeProbes.xlsx'
    lmode_dfA = pd.read_excel(filename, sheet_name='A2', usecols=range(0,14))
    lmode_dfB = pd.read_excel(filename, sheet_name='B2', usecols=range(0,14))
    lmode_dfC = pd.read_excel(filename, sheet_name='C2', usecols=range(0,14))
    net_df = pd.DataFrame()
    net_df['AU R-Rsep (cm)'] = lmode_dfA['rminrsep_U']
    net_df['AD R-Rsep (cm)'] = lmode_dfA['rminrsep_D']
    net_df['BU R-Rsep (cm)'] = lmode_dfB['rminrsep_U']
    net_df['BD R-Rsep (cm)'] = lmode_dfB['rminrsep_D']
    net_df['CU R-Rsep (cm)'] = lmode_dfC['rminrsep_U']
    net_df['CD R-Rsep (cm)'] = lmode_dfC['rminrsep_D']

    # Times 10^4 to go from cm-2 s-1 to m-2 s-1.
    net_df['AU Net Fluence (m-2)'] = lmode_dfA['w_areal_U'] * 10**15 * 10**4
    net_df['AD Net Fluence (m-2)'] = lmode_dfA['w_areal_D'] * 10**15 * 10**4
    net_df['BU Net Fluence (m-2)'] = lmode_dfB['w_areal_U'] * 10**15 * 10**4
    net_df['BD Net Fluence (m-2)'] = lmode_dfB['w_areal_D'] * 10**15 * 10**4
    net_df['CU Net Fluence (m-2)'] = lmode_dfC['w_areal_U'] * 10**15 * 10**4
    net_df['CD Net Fluence (m-2)'] = lmode_dfC['w_areal_D'] * 10**15 * 10**4

    return net_df.dropna(how='all')

def myplot():
     x1   = sput_df['R-Rsep (cm)']
     y1_C = sput_df['Sputt. Fluence from C on 90C10W']
     y1_D = sput_df['Sputt. Fluence from D on 90C10W']
     y1_W = sput_df['Sputt. Fluence from W on 90C10W']
     x1_net = net_df['AD R-Rsep (cm)']
     y1_net = net_df['AD Net Fluence (m-2)']
     x2 = fluxes_df['R-Rsep (cm)']
     y2 = fluxes_df['Te (eV)']
     x3 = fluxes_df['R-Rsep (cm)']
     y3 = fluxes_df['D Flux (m-2 s-1)']
     y4 = fluxes_df['W1+ Flux (m-2 s-1)'] + fluxes_df['W2+ Flux (m-2 s-1)']
     fig = plt.figure()
     fig.suptitle('Comparison of Fluxes on 90% C - 10% W\n0.5% C2+ 0.5% C3+\n0.0005% W1+ 0.0005% W2+')
     ax1 = plt.subplot(411)
     ax2 = plt.subplot(412, sharex=ax1)
     ax3 = plt.subplot(413, sharex=ax1)
     ax4 = plt.subplot(414, sharex=ax1)
     ax1.semilogy(x1, y1_C, '^', alpha=0.4, label='C')
     ax1.semilogy(x1, y1_D, '^', alpha=0.4, label='D')
     ax1.semilogy(x1, y1_W, '^', alpha=0.4, label='W')
     ax1.semilogy(x1, y1_C+y1_D+y1_W, 'C4^', label='C+D+W')
     ax1.semilogy(x1_net, y1_net, 'r-', label='Net')
     ax1.legend()
     ax1.set_ylabel('W Fluence')
     ax2.plot(x2, y2, 'r')
     ax2.set_ylabel('Te (eV)')
     ax3.plot(x3, y3, 'r')
     ax3.set_ylabel('D Flux (m-2 s-1)')
     ax3.set_xlabel('R-Rsep (cm)')
     ax4.semilogy(x3, y4, '>', alpha=0.4, label='Incoming Fluence')
     ax4.semilogy(x3, y1_C+y1_D+y1_W, '>', alpha=0.4,  label='Sputt.')
     ax4.semilogy(x3, y4 - (y1_C+y1_D+y1_W), '>', label='Incoming - Sputt.')
     ax4.semilogy(x1_net, y1_net, 'r-', label='CU Net')
     ax4.set_ylabel('W Fluence')
     ax4.legend()
     ax4.set_xlabel('R-Rsep (cm)')
     ax4.set_ylim([10e15, 10e19])
     #fig.tight_layout()
