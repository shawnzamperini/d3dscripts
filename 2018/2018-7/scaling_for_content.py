from scipy.optimize import curve_fit
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt


filename = '/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets/2018-7/scaling.xlsx'
df = pd.read_excel(filename, sheet_name='Total Content Trend')

# Only use forward H-mode shots.
df.drop([0,1,2,3,9], inplace=True)
# Or use all of them. Comment out whichever.
#df.drop([0,1,2,3], inplace=True)

def fit_func(X1X2, a, b, c):
    x1, x2 = X1X2
    return a * x1**b * x2**c

# Scaling with fg and q95.
if False:

    popt_itf, pcov_itf = curve_fit(fit_func, (df['Greenwald Fraction'], df['q95']), df['ITF per Shot'])
    popt_otf, pcov_otf = curve_fit(fit_func, (df['Greenwald Fraction'], df['q95']), df['OTF per Shot'])

    yfit_itf = fit_func((df['Greenwald Fraction'], df['q95']), *popt_itf)
    yfit_otf = fit_func((df['Greenwald Fraction'], df['q95']), *popt_otf)

    textstr = r'$\mathrm{R_{fit} = a\ *\ f_G^{b}\ *\ Q95^{c}}$' + \
              '\na = {:.2f}'.format(popt_itf[0]) + \
              '\nb = {:.2f}'.format(popt_itf[1]) + \
              '\nc = {:.2f}'.format(popt_itf[2])

    props = dict(alpha=0.25, facecolor='k')
    font = {'fontsize':24, 'weight':'bold'}

    plt.style.use('seaborn')
    font = {'fontsize' : 24, 'weight' : 'bold'}

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(yfit_itf, df['ITF per Shot'], 'k.', label='ITF')
    ax1.set_xlabel('W per Shot (Expected)', font)
    ax1.set_ylabel('W per Shot (Measured)', font)
    ax1.set_xlim([0,0.2])
    ax1.set_ylim([0,0.2])
    ax1.plot(range(0, 5), 'k--')
    ax1.text(0.6, 0.1, textstr, transform=ax1.transAxes, fontsize=18, bbox=props)

    textstr = r'$\mathrm{R_{fit} = a\ *\ f_G^{b}\ *\ Q95^{c}}$' + \
              '\na = {:.2f}'.format(popt_otf[0]) + \
              '\nb = {:.2f}'.format(popt_otf[1]) + \
              '\nc = {:.2f}'.format(popt_otf[2])

    ax2 = fig.add_subplot(212)
    ax2.plot(yfit_otf, df['OTF per Shot'], 'k.', label='OTF')
    ax2.set_xlabel('W per Shot (Expected)')
    ax2.set_ylabel('W per Shot (Measured)')
    ax2.set_xlim([0,0.2])
    ax2.set_ylim([0,0.2])
    ax2.plot(range(0, 5), 'k--')
    ax2.text(0.6, 0.1, textstr, transform=ax2.transAxes, fontsize=18, bbox=props)

    fig.tight_layout()
    fig.show()

# Scaling with Pinj and Ip.
if True:
    # Graph showing PINJ and Ip.

    pinj = df['PINJ (MW)']
    ip   = df['IP (MA)']
    #pinj_norm = np.interp(pinj, (pinj.min(), pinj.max()), (0.0001,1))
    #ip_norm   = np.interp(ip, (ip.min(), ip.max()), (0.0001,1))

    pinj_norm = pinj
    ip_norm = ip

    popt_itf, pcov_itf = curve_fit(fit_func, (pinj_norm, ip_norm), df['ITF per Shot'], bounds=((-np.inf, -2, -2), (np.inf, 2, 2)))
    popt_otf, pcov_otf = curve_fit(fit_func, (pinj_norm, ip_norm), df['OTF per Shot'], bounds=((-np.inf, -2, -2), (np.inf, 2, 2)))

    #popt_itf, pcov_itf = curve_fit(fit_func, (pinj_norm, ip_norm), df['ITF per Shot'], maxfev=5000)
    #popt_otf, pcov_otf = curve_fit(fit_func, (pinj_norm, ip_norm), df['OTF per Shot'], maxfev=5000)

    yfit_itf = fit_func((pinj_norm, ip_norm), *popt_itf)
    yfit_otf = fit_func((pinj_norm, ip_norm), *popt_otf)

    textstr = r'$\mathrm{W_{expected} = a*\ P_{inj}^{b}*I_p^{c}}$' + \
              '\na = {:.2e}'.format(popt_itf[0]) + \
              '\nb = {:.2f}'.format(popt_itf[1]) + \
              '\nc = {:.2f}'.format(popt_itf[2])

    props = dict(alpha=0.25, facecolor='k')
    font = {'fontsize':28, 'weight':'bold'}

    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(yfit_itf, df['ITF per Shot'], 'k.', label='ITF', ms=15)
    ax1.set_xlabel('W per Shot (Expected)', font)
    ax1.set_ylabel('W per Shot (Measured)', font)
    ax1.set_xlim([0,0.2])
    ax1.set_ylim([0,0.2])
    ax1.set_title('ITF', font)
    ax1.tick_params(labelsize=22)
    ax1.plot(range(0, 5), 'k--')
    ax1.text(0.6, 0.1, textstr, transform=ax1.transAxes, fontsize=24, bbox=props)

    textstr = r'$\mathrm{W_{expected} = a*P_{inj}^{b}*I_p^{c}}$' + \
              '\na = {:.2e}'.format(popt_otf[0]) + \
              '\nb = {:.2f}'.format(popt_otf[1]) + \
              '\nc = {:.2f}'.format(popt_otf[2])

    ax2 = fig.add_subplot(122)
    ax2.plot(yfit_otf, df['OTF per Shot'], 'k.', label='OTF', ms=15)
    ax2.set_xlabel('W per Shot (Expected)', font)
    ax2.set_ylabel('W per Shot (Measured)', font)
    ax2.set_xlim([0,0.2])
    ax2.set_ylim([0,0.2])
    ax2.set_title('OTF', font)
    ax2.tick_params(labelsize=22)
    ax2.plot(range(0, 5), 'k--')
    ax2.text(0.6, 0.1, textstr, transform=ax2.transAxes, fontsize=24, bbox=props)

    fig.tight_layout()
    fig.show()

# Scaling with.
if False:
    # Graph showing fg and Ip.

    pinj = df['PINJ (MW)']
    ip   = df['Greenwald Fraction']
    #pinj_norm = np.interp(pinj, (pinj.min(), pinj.max()), (0.0001,1))
    #ip_norm   = np.interp(ip, (ip.min(), ip.max()), (0.0001,1))

    pinj_norm = pinj
    ip_norm = ip

    popt_itf, pcov_itf = curve_fit(fit_func, (pinj_norm, ip_norm), df['ITF per Shot'])
    popt_otf, pcov_otf = curve_fit(fit_func, (pinj_norm, ip_norm), df['OTF per Shot'])

    yfit_itf = fit_func((pinj_norm, ip_norm), *popt_itf)
    yfit_otf = fit_func((pinj_norm, ip_norm), *popt_otf)

    textstr = r'$\mathrm{W_{expected} = a*\ P_{inj}^{b}*f_g^{c}}$' + \
              '\na = {:.2e}'.format(popt_itf[0]) + \
              '\nb = {:.2f}'.format(popt_itf[1]) + \
              '\nc = {:.2f}'.format(popt_itf[2])

    props = dict(alpha=0.25, facecolor='k')
    font = {'fontsize':28, 'weight':'bold'}

    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(yfit_itf, df['ITF per Shot'], 'k.', label='ITF', ms=15)
    ax1.set_xlabel('W per Shot (Expected)', font)
    ax1.set_ylabel('W per Shot (Measured)', font)
    ax1.set_xlim([0,0.2])
    ax1.set_ylim([0,0.2])
    ax1.set_title('ITF', font)
    ax1.tick_params(labelsize=22)
    ax1.plot(range(0, 5), 'k--')
    ax1.text(0.6, 0.1, textstr, transform=ax1.transAxes, fontsize=24, bbox=props)

    textstr = r'$\mathrm{W_{expected} = a*P_{inj}^{b}*f_g^{c}}$' + \
              '\na = {:.2e}'.format(popt_otf[0]) + \
              '\nb = {:.2f}'.format(popt_otf[1]) + \
              '\nc = {:.2f}'.format(popt_otf[2])

    ax2 = fig.add_subplot(122)
    ax2.plot(yfit_otf, df['OTF per Shot'], 'k.', label='OTF', ms=15)
    ax2.set_xlabel('W per Shot (Expected)', font)
    ax2.set_ylabel('W per Shot (Measured)', font)
    ax2.set_xlim([0,0.2])
    ax2.set_ylim([0,0.2])
    ax2.set_title('OTF', font)
    ax2.tick_params(labelsize=22)
    ax2.plot(range(0, 5), 'k--')
    ax2.text(0.6, 0.1, textstr, transform=ax2.transAxes, fontsize=24, bbox=props)

    fig.tight_layout()
    fig.show()

# Scaling with Pinj and Ip, except everything normalized.
if False:
    # Graph showing PINJ and Ip.

    pinj = df['PINJ (MW)']
    ip   = df['IP (MA)']
    pershot_itf = df['ITF per Shot']
    pershot_otf = df['OTF per Shot']
    pinj_norm = np.interp(pinj, (pinj.min(), pinj.max()), (0.0001,1))
    ip_norm   = np.interp(ip, (ip.min(), ip.max()), (0.0001,1))
    itf_norm  = np.interp(pershot_itf, (pershot_itf.min(), pershot_itf.max()), (0.0001, 1))
    otf_norm  = np.interp(pershot_otf, (pershot_otf.min(), pershot_otf.max()), (0.0001, 1))


    popt_itf, pcov_itf = curve_fit(fit_func, (pinj_norm, ip_norm), itf_norm)
    popt_otf, pcov_otf = curve_fit(fit_func, (pinj_norm, ip_norm), otf_norm)

    yfit_itf = fit_func((pinj_norm, ip_norm), *popt_itf)
    yfit_otf = fit_func((pinj_norm, ip_norm), *popt_otf)

    textstr = r'$\mathrm{W_{expected} = a*\ P_{inj}^{b}*I_p^{c}}$' + \
              '\na = {:.2e}'.format(popt_itf[0]) + \
              '\nb = {:.2f}'.format(popt_itf[1]) + \
              '\nc = {:.2f}'.format(popt_itf[2])

    props = dict(alpha=0.25, facecolor='k')
    font = {'fontsize':28, 'weight':'bold'}

    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(yfit_itf, itf_norm, 'k.', label='ITF', ms=15)
    ax1.set_xlabel('W per Shot (Expected)', font)
    ax1.set_ylabel('W per Shot (Measured)', font)
    ax1.set_xlim([0,1])
    ax1.set_ylim([0,1])
    ax1.set_title('ITF', font)
    ax1.tick_params(labelsize=22)
    ax1.plot(range(0, 5), 'k--')
    ax1.text(0.6, 0.1, textstr, transform=ax1.transAxes, fontsize=24, bbox=props)

    textstr = r'$\mathrm{W_{expected} = a*P_{inj}^{b}*I_p^{c}}$' + \
              '\na = {:.2e}'.format(popt_otf[0]) + \
              '\nb = {:.2f}'.format(popt_otf[1]) + \
              '\nc = {:.2f}'.format(popt_otf[2])

    ax2 = fig.add_subplot(122)
    ax2.plot(yfit_otf, otf_norm, 'k.', label='OTF', ms=15)
    ax2.set_xlabel('W per Shot (Expected)', font)
    ax2.set_ylabel('W per Shot (Measured)', font)
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    ax2.set_title('OTF', font)
    ax2.tick_params(labelsize=22)
    ax2.plot(range(0, 5), 'k--')
    ax2.text(0.6, 0.1, textstr, transform=ax2.transAxes, fontsize=24, bbox=props)

    fig.tight_layout()
    fig.show()

# Scaling with Pinj and q95.
if False:

    pinj = df['PINJ (MW)']
    q95   = df['q95']
    #pinj_norm = np.interp(pinj, (pinj.min(), pinj.max()), (0.0001,1))
    #ip_norm   = np.interp(ip, (ip.min(), ip.max()), (0.0001,1))

    popt_itf, pcov_itf = curve_fit(fit_func, (pinj, q95), df['ITF per Shot'], maxfev=5000)
    popt_otf, pcov_otf = curve_fit(fit_func, (pinj, q95), df['OTF per Shot'], maxfev=5000)

    yfit_itf = fit_func((pinj, q95), *popt_itf)
    yfit_otf = fit_func((pinj, q95), *popt_otf)

    textstr = r'$\mathrm{W_{expected} = a*\ P_{inj}^{b}*q_{95}^{c}}$' + \
              '\na = {:.2e}'.format(popt_itf[0]) + \
              '\nb = {:.2f}'.format(popt_itf[1]) + \
              '\nc = {:.2f}'.format(popt_itf[2])

    props = dict(alpha=0.25, facecolor='k')
    font = {'fontsize':28, 'weight':'bold'}

    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(yfit_itf, df['ITF per Shot'], 'k.', label='ITF', ms=15)
    ax1.set_xlabel('W per Shot (Expected)', font)
    ax1.set_ylabel('W per Shot (Measured)', font)
    ax1.set_xlim([0,0.2])
    ax1.set_ylim([0,0.2])
    ax1.set_title('ITF', font)
    ax1.tick_params(labelsize=22)
    ax1.plot(range(0, 5), 'k--')
    ax1.text(0.6, 0.05, textstr, transform=ax1.transAxes, fontsize=24, bbox=props)

    textstr = r'$\mathrm{W_{expected} = a*P_{inj}^{b}*q_{95}^{c}}$' + \
              '\na = {:.2e}'.format(popt_otf[0]) + \
              '\nb = {:.2f}'.format(popt_otf[1]) + \
              '\nc = {:.2f}'.format(popt_otf[2])

    ax2 = fig.add_subplot(122)
    ax2.plot(yfit_otf, df['OTF per Shot'], 'k.', label='OTF', ms=15)
    ax2.set_xlabel('W per Shot (Expected)', font)
    ax2.set_ylabel('W per Shot (Measured)', font)
    ax2.set_xlim([0,0.2])
    ax2.set_ylim([0,0.2])
    ax2.set_title('OTF', font)
    ax2.tick_params(labelsize=22)
    ax2.plot(range(0, 5), 'k--')
    ax2.text(0.6, 0.05, textstr, transform=ax2.transAxes, fontsize=24, bbox=props)

    fig.tight_layout()
    fig.show()

if False:

    def fit_func(X1X2X3, a, b, c, d):
        x1, x2, x3 = X1X2X3
        return a * x1**b * x2**c * x3**d


    pinj = df['PINJ (MW)']
    q95   = df['q95']
    lamb = df['Density Fall Off (cm)']
    #pinj_norm = np.interp(pinj, (pinj.min(), pinj.max()), (0.0001,1))
    #ip_norm   = np.interp(ip, (ip.min(), ip.max()), (0.0001,1))

    popt_itf, pcov_itf = curve_fit(fit_func, (pinj, q95, lamb), df['ITF per Shot'], maxfev=5000)
    popt_otf, pcov_otf = curve_fit(fit_func, (pinj, q95, lamb), df['OTF per Shot'], maxfev=5000)

    yfit_itf = fit_func((pinj, q95, lamb), *popt_itf)
    yfit_otf = fit_func((pinj, q95, lamb), *popt_otf)

    textstr = r'$\mathrm{W_{expected} = a\ P_{inj}^{b}\ q_{95}^{c}\ \lambda_{ne}^{d}}$' + \
              '\na = {:.2e}'.format(popt_itf[0]) + \
              '\nb = {:.2f}'.format(popt_itf[1]) + \
              '\nc = {:.2f}'.format(popt_itf[2]) + \
              '\nd = {:.2f}'.format(popt_itf[3])

    props = dict(alpha=0.25, facecolor='k')
    font = {'fontsize':28, 'weight':'bold'}

    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(yfit_itf, df['ITF per Shot'], 'k.', label='ITF', ms=15)
    ax1.set_xlabel('W per Shot (Expected)', font)
    ax1.set_ylabel('W per Shot (Measured)', font)
    ax1.set_xlim([0,0.2])
    ax1.set_ylim([0,0.2])
    ax1.set_title('ITF', font)
    ax1.tick_params(labelsize=22)
    ax1.plot(range(0, 5), 'k--')
    ax1.text(0.6, 0.05, textstr, transform=ax1.transAxes, fontsize=20, bbox=props)

    textstr = r'$\mathrm{W_{expected} = a\ P_{inj}^{b}\ q_{95}^{c}\ \lambda_{ne}^{d}}$' + \
              '\na = {:.2e}'.format(popt_otf[0]) + \
              '\nb = {:.2f}'.format(popt_otf[1]) + \
              '\nc = {:.2f}'.format(popt_otf[2]) + \
              '\nd = {:.2f}'.format(popt_otf[3])

    ax2 = fig.add_subplot(122)
    ax2.plot(yfit_otf, df['OTF per Shot'], 'k.', label='OTF', ms=15)
    ax2.set_xlabel('W per Shot (Expected)', font)
    ax2.set_ylabel('W per Shot (Measured)', font)
    ax2.set_xlim([0,0.2])
    ax2.set_ylim([0,0.2])
    ax2.set_title('OTF', font)
    ax2.tick_params(labelsize=22)
    ax2.plot(range(0, 5), 'k--')
    ax2.text(0.6, 0.05, textstr, transform=ax2.transAxes, fontsize=20, bbox=props)

    fig.tight_layout()
    fig.show()
