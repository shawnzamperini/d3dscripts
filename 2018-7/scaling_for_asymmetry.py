from scipy.optimize      import curve_fit
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt


# Load everything into a Dataframe (the "Already RBS'd" part).
filename = '/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets/2018-7/scaling.xlsx'
df = pd.read_excel(filename, sheet_name='Scaling Law 1', skiprows=12)

# Drop everything after, except the three unique ones.
df.drop([9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24], inplace=True)

def scale_fg_lamb():

    ratios = np.array(df['Ratio'].values, dtype=np.float64)
    fgs = np.array(df['Greenwald Fraction'].values, dtype=np.float64)
    lambs = np.array(df['Density Fall Off (cm)'].values, dtype=np.float64)

    # Rescale the independent variables to 0-1.
    #ratios_scale = np.interp(ratios, (ratios.min(), ratios.min()), (0,1))
    fgs_scale = np.interp(fgs, (fgs.min(), fgs.max()), (0.0001,1))
    lambs_scale = np.interp(lambs, (lambs.min(), lambs.max()), (0.0001,1))

    def fit_func(X1X2, a, b, c):
        x1, x2 = X1X2
        return a * x1**b * x2**c

    popt, pcov = curve_fit(fit_func, (fgs_scale, lambs_scale), ratios)
    fit_vals = fit_func((fgs_scale, lambs_scale), *popt)

    textstr = r'$\mathrm{R_{scaling} = A\ *\ f_g^{B}\ *\ \lambda_{ne}^{C}}$' + \
              '\nA = {:.2f}'.format(popt[0]) + \
              '\nB = {:.2f}'.format(popt[1]) + \
              '\nC = {:.2f}'.format(popt[2])

    # Some properties for the plots.
    props = dict(alpha=0.25, facecolor='k')
    font = {'fontsize':24, 'weight':'bold'}
    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(fit_vals, ratios, 'k.', ms=10)
    ax1.set_xlabel(r'$\mathrm{R_{scaling}}$', font)
    ax1.set_ylabel(r'$\mathrm{R_{measured}}$', font)
    ax1.plot(range(0,5), 'k--')
    ax1.set_xlim([0,2])
    ax1.set_ylim([0,2])
    ax1.set_title('ITF/OTF Total Content Ratio', font)
    plt.tick_params(axis='both', which='major', labelsize=18)
    ax1.text(0.6, 0.1, textstr, transform=ax1.transAxes, fontsize=18, bbox=props)
    fig.tight_layout()
    fig.show()

    return fit_vals, ratios
