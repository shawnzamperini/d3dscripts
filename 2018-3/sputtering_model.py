import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import adf11_util_3 as adf11
import warnings
import calc_rsep as rsep
from scipy.signal import savgol_filter
from scipy.interpolate import interpolate



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
    #E0 = [float(E) for E in E0]

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
    # Weird warning only involving np.floats, but doesn't seem to cause problems so ignore it.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
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

def yields(Emin=0.001, Emax=1200):
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

def get_fluxes(#carbon_frac=  np.array([0.001,  0.001,  0.001,  0.001,  0.001,  0.001]) / 6.0,
               #tungsten_frac=np.array([0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001]) / 6.0):
               carbon_frac, tungsten_frac):
    """
    Calculates the fluxes of each ion charge state. The elements in the lists
    passed in each correspond to a charge state (i.e. C1+, C2+, ..., C6+).

    Default is [...] / 6.0 because if we assume 1% of the flux is C, then if we
    want to evenly distribute it across the charge states, divide by the 6 states.
    A better analysis would come from a model like DIVIMP or something.
    """

    excel_file = "/home/shawn/d3dscripts/Data/LP_with_fit.xlsx"

    print("Tungsten fraction: " + str(tungsten_frac))
    # Use plunge 2 of 167195 since it gets the closest to the separatrix.
    #df_195_2 = pd.read_excel(excel_file,
    #                         sheet_name="LP Data",
    #                         skiprows=[0, 1],
    #                         names=["Time (ms)", "R (cm)", "ne (e18 m-3)", "Te (eV)",
    #                                "Vfl (V)", "Vp (V)", "R-Rsep (cm)"],
    #                         usecols=[57, 58, 59, 60, 61, 62, 63])

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

    # Try using all the data.
    if True:
        lp_df = pd.read_excel(excel_file,
                              sheet_name="Fit Te's",
                              #skiprows=[0],
                              names=['R (cm)', 'Te (eV)', 'ne (e18 m-3)'],
                              usecols=[0, 1, 3])

    # Want the flow at probe face (= flow at sheath edge). It's 0.5ne*cs
    m_deut = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
    cs_series = (2 * lp_df["Te (eV)"] / m_deut) ** (1/2)
    flux_series = 0.5 * lp_df['ne (e18 m-3)'] * cs_series * 10**(18)
    lp_df['D Flux (m-2 s-1)'] = flux_series

    # Estimate of the carbon and tungsten flux as fractions of the deuterium flux.
    for charge_state in range(0, 6):
        lp_df['C' + str(charge_state+1) +  '+ Flux (m-2 s-1)'] = \
            lp_df['D Flux (m-2 s-1)'] * carbon_frac[charge_state]
        lp_df['W' + str(charge_state+1) + '+ Flux (m-2 s-1)'] = \
            lp_df['D Flux (m-2 s-1)'] * tungsten_frac[charge_state]

    return lp_df

def get_sput_flux(fluxes_df):
    y_df = yields()
    #fluxes_df = get_fluxes()
    sput_df = pd.DataFrame()
    sput_df["R (cm)"]      = fluxes_df["R (cm)"]
    #sput_df["R-Rsep (cm)"] = fluxes_df["R-Rsep (cm)"]

    # Calculate the sputtered flux from each charge state of the ions.
    for ion in ["Carbon", "Tungsten", "Deuterium"]:
        # As of now only considering charge states 1-6.
        for charge_state in range(1, 7):
            if ion == "Carbon":
                fluxes = fluxes_df["C" +str(charge_state) + "+ Flux (m-2 s-1)"].values
                Tes    = fluxes_df["Te (eV)"].values
                sput_flux = np.array([])
                for flux, Te in zip(fluxes, Tes):
                    tmp_y    = get_yield(y_df, ion.lower(), charge_state, Te, verbose=False)
                    tmp_sput = tmp_y * flux
                    sput_flux = np.append(sput_flux, tmp_sput)
                sput_df['Sputt. Flux from C' + str(charge_state) + '+'] = sput_flux
            elif ion == "Tungsten":
                fluxes = fluxes_df["W" + str(charge_state) + "+ Flux (m-2 s-1)"].values
                Tes    = fluxes_df["Te (eV)"].values
                sput_flux = np.array([])
                for flux, Te in zip(fluxes, Tes):
                    tmp_y    = get_yield(y_df, ion.lower(), charge_state, Te, verbose=False)
                    tmp_sput = tmp_y * flux
                    sput_flux = np.append(sput_flux, tmp_sput)
                sput_df['Sputt. Flux from W' + str(charge_state) + '+'] = sput_flux
            elif ion == "Deuterium":
                if charge_state > 1:
                    pass
                else:
                    fluxes = fluxes_df["D Flux (m-2 s-1)"].values
                    Tes    = fluxes_df["Te (eV)"].values
                    sput_flux = np.array([])
                    for flux, Te in zip(fluxes, Tes):
                        tmp_y    = get_yield(y_df, ion.lower(), charge_state, Te, verbose=False)
                        tmp_sput = tmp_y * flux
                        sput_flux = np.append(sput_flux, tmp_sput)
                    sput_df['Sputt. Flux from D'] = sput_flux

    # Get the total sputtered flux from each ion by adding up each contribution.
    sput_df["Sputt. Flux from C"] = sput_df["Sputt. Flux from C1+"].values + \
                                    sput_df["Sputt. Flux from C2+"].values + \
                                    sput_df["Sputt. Flux from C3+"].values + \
                                    sput_df["Sputt. Flux from C4+"].values + \
                                    sput_df["Sputt. Flux from C5+"].values + \
                                    sput_df["Sputt. Flux from C6+"].values
    sput_df["Sputt. Flux from W"] = sput_df["Sputt. Flux from W1+"].values + \
                                    sput_df["Sputt. Flux from W2+"].values + \
                                    sput_df["Sputt. Flux from W3+"].values + \
                                    sput_df["Sputt. Flux from W4+"].values + \
                                    sput_df["Sputt. Flux from W5+"].values + \
                                    sput_df["Sputt. Flux from W6+"].values

    return sput_df

def plot_sput_flux(fluxes_df):
    sput_df = get_sput_flux(fluxes_df)
    x = sput_df["R-Rsep (cm)"]
    yC = sput_df["Sputt. Flux from C"].values
    yW = sput_df["Sputt. Flux from W"].values
    yD = sput_df["Sputt. Flux from D"].values

    plt.semilogy(x, yC, label='carbon 1%')
    plt.semilogy(x, yW, label='tungsten 0.1%')
    plt.semilogy(x, yD, label='deuterium')
    plt.legend()
    plt.xlabel("R-Rsep (cm)")
    plt.ylabel(r"$\mathrm{Sputt. Flux (m^{-2}\ s^{-1)}}$")
    plt.show()

def calc_iz_frac(lp_df, filename='/home/shawn/d3dscripts/Data/adf11/scd50/scd50_w.dat'):
    # Probe widths in m
    a_size = np.float64(3.0 / 100.0)
    b_size = np.float64(1.0 / 100.0)
    c_size = np.float64(0.5 / 100.0)
    # Mass of tungsten in eV s^2 m^-2.
    mass_w = 183.84 * 931.49 * 10**6.0 / ((3*10**8.0)**2.0)

    print("Loading ADAS file: " + filename)
    a = adf11.adf11(filename, debug=False)
    a.interpolate()
    # lambda_iz = speed_w * (ne * sigmanu_bar)^-1
    # Let's get everything separately and in the right units.
    # In eV.
    Te = lp_df['Te (eV)']

    # In e18 m-3.
    ne = lp_df['ne (e18 m-3)']
    # In m-3
    ne = ne * 10**18

    # In m/s.
    v_w = np.sqrt(3 * Te / mass_w)

    # In cm3/s. 0th ionization state since we assume anything > 1 will return.
    sigmanu = np.array([])
    for ne_tmp, te_tmp in zip(ne, Te):
        sigmanu_tmp = a.EvalInterpolation(ne_tmp, te_tmp, 0)
        # In m3/s.
        sigmanu_tmp = sigmanu_tmp * 10**(-6)
        sigmanu = np.append(sigmanu, sigmanu_tmp)

    # Finally, lambda_iz in m.
    lambda_iz = v_w * (ne * sigmanu)**(-1)

    # Put into lp_df.
    lp_df['Lambda_iz (m)'] = lambda_iz

    # Calculate fraction that ionizes within each probe tube.
    lp_df['A Frac. Ion.'] = 1.0 - np.exp(- a_size / lp_df['Lambda_iz (m)'])
    lp_df['B Frac. Ion.'] = 1.0 - np.exp(- b_size / lp_df['Lambda_iz (m)'])
    lp_df['C Frac. Ion.'] = 1.0 - np.exp(- c_size / lp_df['Lambda_iz (m)'])

    return lp_df

def calc_lost_flux(lp_df, sput_df):
    # 1 - frac_ion bc that which does not ionize is lost.
    a_lost_c = sput_df['Sputt. Flux from C'] * (1.0 - lp_df['A Frac. Ion.'])
    b_lost_c = sput_df['Sputt. Flux from C'] * (1.0 - lp_df['B Frac. Ion.'])
    c_lost_c = sput_df['Sputt. Flux from C'] * (1.0 - lp_df['C Frac. Ion.'])
    a_lost_w = sput_df['Sputt. Flux from W'] * (1.0 - lp_df['A Frac. Ion.'])
    b_lost_w = sput_df['Sputt. Flux from W'] * (1.0 - lp_df['B Frac. Ion.'])
    c_lost_w = sput_df['Sputt. Flux from W'] * (1.0 - lp_df['C Frac. Ion.'])

    lost_df = pd.DataFrame()
    lost_df['A Lost Flux from C'] = a_lost_c
    lost_df['B Lost Flux from C'] = b_lost_c
    lost_df['C Lost Flux from C'] = c_lost_c
    lost_df['A Lost Flux from W'] = a_lost_w
    lost_df['B Lost Flux from W'] = b_lost_w
    lost_df['C Lost Flux from W'] = c_lost_w

    # Odds and ends to pass along.
    lost_df['R (cm)']       = lp_df['R (cm)']
    lost_df['Te (eV)']      = lp_df['Te (eV)']
    lost_df['ne (e18 m-3)'] = lp_df['ne (e18 m-3)']
    return lost_df

def calc_net_flux():
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
    net_df['AU Net Flux (m-2 s-1)'] = lmode_dfA['w_areal_U'] * 10**15 * 10**4 / exposed_time
    net_df['AD Net Flux (m-2 s-1)'] = lmode_dfA['w_areal_D'] * 10**15 * 10**4 / exposed_time
    net_df['BU Net Flux (m-2 s-1)'] = lmode_dfB['w_areal_U'] * 10**15 * 10**4 / exposed_time
    net_df['BD Net Flux (m-2 s-1)'] = lmode_dfB['w_areal_D'] * 10**15 * 10**4 / exposed_time
    net_df['CU Net Flux (m-2 s-1)'] = lmode_dfC['w_areal_U'] * 10**15 * 10**4 / exposed_time
    net_df['CD Net Flux (m-2 s-1)'] = lmode_dfC['w_areal_D'] * 10**15 * 10**4 / exposed_time
    net_df['AU W Areal (1e19 m-2)'] = lmode_dfA['w_areal_U']
    net_df['AD W Areal (1e19 m-2)'] = lmode_dfA['w_areal_D']
    net_df['BU W Areal (1e19 m-2)'] = lmode_dfB['w_areal_U']
    net_df['BD W Areal (1e19 m-2)'] = lmode_dfB['w_areal_D']
    net_df['CU W Areal (1e19 m-2)'] = lmode_dfC['w_areal_U']
    net_df['CD W Areal (1e19 m-2)'] = lmode_dfC['w_areal_D']


    return net_df

def plot_net_lost(net_df, lost_df, probe='AD'):
    x1 = net_df[probe + " R-Rsep (cm)"]
    y1 = net_df[probe + " Net Flux (m-2 s-1)"]
    x2 = lost_df["A R-Rsep (cm)"]
    y2_c = lost_df['A Lost Flux from C']
    y2_w = lost_df['A Lost Flux from W']

    plt.rcParams.update({'font.size': 34})
    plt.plot(x1, y1, 'b', label='Net Flux to Probe')
    plt.plot(x2, y2_c, 'r', label='Lost Flux off Probe (C)')
    plt.plot(x2, y2_w, 'r--', label='Lost Flux off Probe (W)')
    plt.xlabel('R-Rsep (cm)')
    plt.ylabel(r"$\mathrm{Flux (m^{-2}\ s^{-1)}}$")
    plt.title(probe + ' Fluxes due to Sputtering')
    plt.legend()
    plt.show()

def plot_areal_lost(net_df, lost_df, probe):
    x1        = lost_df[probe[0] + ' R-Rsep (cm)']
    y1        = lost_df[probe[0] + ' Lost Flux from C']
    y1_smooth = savgol_filter(y1, 21, 2)
    x2        = net_df[probe + " R-Rsep (cm)"]
    y2        = net_df[probe + ' W Areal (1e19 m-2)']

    # Create a common basis.
    f_y1 = interpolate.interp1d(x1, y1_smooth)
    f_y2 = interpolate.interp1d(x2, y2)
    common_x = np.linspace(np.max([x1.min(), x2.min()]), np.min([x1.max(), x2.max()]), 100)
    new_y1 = f_y1(common_x)
    new_y2 = f_y2(common_x)

    # Threshold where sputtered flux is greater than 1% of the net flux.
    thresh_low = np.argmin(new_y1 > 0.01 * new_y2 * 10**19)
    thresh_high = np.argmin(new_y1 > 0.1 * new_y2 * 10**19)
    #print("Threshold met at {0}".format(common_x[thresh]))
    plt.rcParams.update({'font.size': 42})
    plt.rcParams.update({'figure.autolayout': True})
    #plt.legend(prop={"size":26})
    fig1, ax1 = plt.subplots()
    ax1.plot(common_x, new_y1, 'r')
    ax1.set_xlabel('R-Rsep (cm)')
    ax1.set_ylabel('Sputtered Flux (m-2 s-1)', color='r')
    ax1.tick_params('y', colors='r')
    ax2 = ax1.twinx()
    ax2.plot(common_x, new_y2, 'b')
    ax2.set_ylabel('W Areal Density (1e19 m-2)', color='b')
    ax2.tick_params('y', colors='b')
    #ax1.axvline(x=common_x[thresh])
    ax1.axvspan(common_x[thresh_low], common_x[thresh_high], alpha = 0.5, color='red')
    plt.show()

def plot_smoothing(net_df, lost_df, lp_df, probe):
    x1 = net_df[probe + " R-Rsep (cm)"]
    y1 = net_df[probe + ' W Areal (1e19 m-2)']

    # Data we'd like to smooth.
    x2 = lost_df["A R-Rsep (cm)"]
    y2 = lost_df['A Lost Flux from C']
    y2_smoothed = savgol_filter(y2, 31, 2, mode='nearest')
    y3 = lp_df[probe[0] + ' Frac. Ion.']

    # Now let's plot it.
    plt.rcParams.update({'font.size': 34})
    fig1, ax1 = plt.subplots()
    ax1.plot(x1, y1, 'b', label='Areal Density')
    ax1.set_xlabel('R-Rsep (cm)')
    ax1.set_ylabel('W Areal Density (1e19 m-2)', color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.plot(x2, y2_smoothed, 'r', label='Lost Flux')
    ax2.set_ylabel('Lost Flux (m-2 s-1)', color='r')
    ax2.tick_params('y', colors='r')
    plt.title(probe + ' Fluxes due to Sputtering')
    plt.legend()
    plt.show()

    plt.plot(x2, y3, '--')
    plt.show()

def plot_sputt_and_areal(sput, net, probe):
    x1        = sput[probe[0] + ' R-Rsep (cm)']
    y1        = sput['Sputt. Flux from C']
    y1_smooth = savgol_filter(y1, 11, 3)
    y1_smoother = savgol_filter(y1, 91, 3)
    x2        = net[probe + ' R-Rsep (cm)']
    y2        = net[probe + ' W Areal (1e19 m-2)']

    # Create a common basis.
    np.warnings.filterwarnings('ignore')
    f_y1 = interpolate.interp1d(x1, y1_smooth)
    f_y1_smoother = interpolate.interp1d(x1, y1_smoother)
    f_y2 = interpolate.interp1d(x2, y2)
    common_x = np.linspace(np.max([x1.min(), x2.min()]), np.min([x1.max(), x2.max()]), 100)
    new_y1 = f_y1(common_x)
    new_y1_smoother = f_y1_smoother(common_x)
    new_y2 = f_y2(common_x)

    # Threshold where sputtered flux is greater than 1% of the net flux.
    thresh_low = np.argmin(new_y1 > 0.01 * new_y2 * 10**19)
    thresh_high = np.argmin(new_y1 > 0.1 * new_y2 * 10**19)
    #print("Threshold met at {0}".format(common_x[thresh]))
    plt.rcParams.update({'font.size': 42})
    plt.rcParams.update({'figure.autolayout': True})
    #plt.legend(prop={"size":26})
    fig1, ax1 = plt.subplots()
    ax1.plot(common_x, new_y1, 'r.')
    ax1.plot(common_x, new_y1_smoother, 'r')
    ax1.set_xlabel('R-Rsep (cm)')
    ax1.set_ylabel('Sputtered Flux (m-2 s-1)', color='r')
    ax1.tick_params('y', colors='r')
    ax2 = ax1.twinx()
    ax2.plot(common_x, new_y2, 'b')
    ax2.set_ylabel('W Areal Density (1e19 m-2)', color='b')
    ax2.tick_params('y', colors='b')
    #ax1.axvline(x=common_x[thresh])
    ax1.axvspan(common_x[thresh_low], common_x[thresh_high], alpha = 0.5, color='red')
    plt.show()

def run_script(plotting=None, get_rsep=True, time_start=1500, time_end=4500, time_step=500):

    #                         C1+     C2+     ...                     C6+
    carbon_frac   = np.array([0.000,  0.005,  0.005,  0.000,  0.000,  0.000])
    tungsten_frac = np.array([0.000, 0.000, 0.000, 0.000, 0.000, 0.000])
    print()
    print("Calculating fluxes of each charge state as fractions of D flux...")
    print("              1+       2+       3+       4+       5+       6+")
    print("C Fraction    {0:3.5f}  {1:3.5f}  {2:3.5f}  {3:3.5f}  {4:3.5f}  {5:3.5f}"
          .format(carbon_frac[0], carbon_frac[1], carbon_frac[2], carbon_frac[3],
                  carbon_frac[4], carbon_frac[5]))
    print("W Fraction    {0:3.5f}  {1:3.5f}  {2:3.5f}  {3:3.5f}  {4:3.5f}  {5:3.5f}"
          .format(tungsten_frac[0], tungsten_frac[1], tungsten_frac[2], tungsten_frac[3],
                  tungsten_frac[4], tungsten_frac[5]))
    print()
    lp_df = get_fluxes(carbon_frac, tungsten_frac)

    print("Calculating sputtered flux...")
    sput_df = get_sput_flux(lp_df)
    print("Calculating lost flux...")
    lp_df = calc_iz_frac(lp_df)
    lost_df = calc_lost_flux(lp_df, sput_df)

    print("Calculating average Rsep for 167192 in range " + str(time_start) + "-" + str(time_end) + " ms")
    if get_rsep == True:
        rsep_dict = rsep.return_avg(167192, time_start, time_end, time_step)
        # Note that LP data is ~ Z location of A probe only. B and C are approximations without
        # mapping to psin and such.
        lost_df['A R-Rsep (cm)'] = lost_df['R (cm)'] - rsep_dict['Average Rsep at A Probe'] * 100.0
        lost_df['B R-Rsep (cm)'] = lost_df['R (cm)'] - rsep_dict['Average Rsep at B Probe'] * 100.0
        lost_df['C R-Rsep (cm)'] = lost_df['R (cm)'] - rsep_dict['Average Rsep at C Probe'] * 100.0
        sput_df['A R-Rsep (cm)'] = sput_df['R (cm)'] - rsep_dict['Average Rsep at A Probe'] * 100.0
        sput_df['B R-Rsep (cm)'] = sput_df['R (cm)'] - rsep_dict['Average Rsep at B Probe'] * 100.0
        sput_df['C R-Rsep (cm)'] = sput_df['R (cm)'] - rsep_dict['Average Rsep at C Probe'] * 100.0

    print("Calculating net flux...")
    net_df = calc_net_flux()
    print("Plotting...")
    if plotting == 1:
        plot_net_lost(net_df, lost_df, probe='AD')
        plot_net_lost(net_df, lost_df, probe='AU')
        plot_net_lost(net_df, lost_df, probe='BD')
        plot_net_lost(net_df, lost_df, probe='BU')
        plot_net_lost(net_df, lost_df, probe='CD')
        plot_net_lost(net_df, lost_df, probe='CU')
    elif plotting == 2:
        plot_areal_lost(net_df, lost_df, probe='AD')
        plot_areal_lost(net_df, lost_df, probe='AU')
        plot_areal_lost(net_df, lost_df, probe='BD')
        plot_areal_lost(net_df, lost_df, probe='BU')
        plot_areal_lost(net_df, lost_df, probe='CD')
        plot_areal_lost(net_df, lost_df, probe='CU')

    print("\nDone.")


    return lp_df, sput_df, lost_df, net_df
