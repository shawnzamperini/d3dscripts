import pandas as pd
import numpy  as np
import pretty_plots as pp
import matplotlib.pyplot as plt


top_path = '~/Drive/School/Tennessee/Research/Collector Probe Excel Sheets/'

# Probes with the respective shot that has MAFOT data.
probes = {'A2' :167196,
          'A3' :167227,
          'A4' :167229,
          'A7' :167237,
          'A8' :167247,
          'A11':167266,
          'A12':167268,
          'A15':167277,
          'A17':167279,
          'A18':167320,
          'A19':167321,
          'A20':167322,
          'A21':167353,
          'A28':167405,
          'A33':167530,
          'A34':167534,
          'A35':167536}

def run(probe, psin_offset=0.0, identify_psin=True):

    # First load the probe Excel sheet.
    rbs_df = pd.read_excel(top_path + probe + '.xlsx')
    shot = probes[probe]

    # Get the psin and W data, grabbing the correct side.
    if probe in ['A2', 'A3', 'A4', 'A7', 'A8']:
        itf_side = 'D'
        otf_side = 'U'
    else:
        itf_side = 'U'
        otf_side = 'D'

    psin_itf = rbs_df['Psin ' + itf_side]
    psin_otf = rbs_df['Psin ' + otf_side]
    w_itf    = rbs_df['W Areal Density ' + itf_side + ' (1e15 W/cm2)']
    w_otf    = rbs_df['W Areal Density ' + otf_side + ' (1e15 W/cm2)']

    # Then grab the MAFOT data. Need to make sure I keep these in the same format.
    mafot_path = top_path + 'Connection Lengths/' + str(shot) + '/' + str(shot) + '.xlsx'
    mafot_itf_df = pd.read_excel(mafot_path, sheet_name='MAFOT ITF', skiprows=2)
    mafot_otf_df = pd.read_excel(mafot_path, sheet_name='MAFOT OTF', skiprows=2)

    # Clean up the data a bit.
    psin_itf = psin_itf[~np.isnan(psin_itf)]
    psin_otf = psin_otf[~np.isnan(psin_otf)]
    w_itf = w_itf[~np.isnan(w_itf)]
    w_otf = w_otf[~np.isnan(w_otf)]

    # Find closest index in mafot data to this psin.
    psin_mafot_itf = mafot_itf_df['psi'].values
    psin_mafot_otf = mafot_otf_df['psi'].values
    idx_itf = [(np.abs(psin_mafot_itf - psin + psin_offset)).argmin() for psin in psin_itf]
    idx_otf = [(np.abs(psin_mafot_otf - psin + psin_offset)).argmin() for psin in psin_otf]

    # Grab the connection lengths, convert from km to cm.
    itf_L = mafot_itf_df['Connection Length (km)'].values[idx_itf] * 100000
    otf_L = mafot_otf_df['Connection Length (km)'].values[idx_otf] * 100000

    # We may need help identifying what the psin_offset should be. Let's plot
    # it out on a log plot and see if it may be obvious.
    if identify_psin:
        fig1, ax1 = plt.subplots()
        ax1.set_xlabel('Psin')
        ax1.set_ylabel('W Areal Density', color='r')
        ax1.semilogy(psin_itf, w_itf, '.', color='r')
        ax1.tick_params(axis='y', labelcolor='r')
        ax1.set_xlim([1.1, 1.4])
        #ax1.set_ylim([5e2, None])
        ax2 = ax1.twinx()
        ax2.set_ylabel('Lconn (cm)')
        ax2.semilogy(psin_mafot_itf+psin_offset, mafot_itf_df['Connection Length (km)'] * 100000, '-', color='k')
        ax2.set_title('ITF')
        fig1.tight_layout()
        fig1.show()

        fig2, ax3 = plt.subplots()
        ax3.set_xlabel('Psin')
        ax3.set_ylabel('W Areal Density', color='r')
        ax3.semilogy(psin_otf, w_otf, '.', color='r')
        ax3.tick_params(axis='y', labelcolor='r')
        ax3.set_xlim([1.1, 1.4])
        #ax3.set_ylim([5e2, None])
        ax4 = ax3.twinx()
        ax4.set_ylabel('Lconn (cm)')
        ax4.semilogy(psin_mafot_otf+psin_offset, mafot_otf_df['Connection Length (km)'] * 100000, '-', color='k')
        ax4.set_title('OTF')
        #ax4.set_ylim([1e1, 1e3])
        fig2.tight_layout()
        fig2.show()

    # Calculate the W per cm.
    w_per_cm_itf = w_itf / itf_L
    w_per_cm_otf = w_otf / otf_L

    # Plot it out.
    fig = pp.pplot(psin_itf, w_per_cm_itf, color=8, label='ITF', fmt='-')
    fig = pp.pplot(psin_otf, w_per_cm_otf, color=6, label='OTF', fmt='-', xlabel='Psin', ylabel='W per Lconn', fig=fig)
