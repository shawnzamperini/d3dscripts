import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt
from scipy import interpolate


# First load the Excel sheet data into pandas DataFrame
df = pd.read_excel('/home/shawn/d3dscripts/Data/LPData.xlsx',
                         sheet_name='Lengths',
                         skiprows=3,
                         usecols=[4,5,6,11,17,18,19,23,24,25,37],
                         names=['A R-Rsep omp (cm)', 'B R-Rsep omp (cm)',
                                 'C R-Rsep omp (cm)', 'Conn R-Rsep omp (cm)',
                                 'A ODF Conn (m)'   , 'B ODF Conn (m)',
                                 'C ODF Conn (m)'   , 'A IDF Conn (m)',
                                 'B IDF Conn (m)'   , 'C IDF Conn (m)',
                                 'Sound Speed (m/s)'])

# Interpolations of conn(r-rsep omp). The isnan is to filter out the nan's (the
# ~ is a negation, so it indexes True with all the not nans).
tmp_x = df['Conn R-Rsep omp (cm)']
tmp_y = df['A ODF Conn (m)']
f_conn_odf_a = interpolate.interp1d(tmp_x[~np.isnan(tmp_x)], tmp_y[~np.isnan(tmp_y)])
tmp_x = df['Conn R-Rsep omp (cm)']
tmp_y = df['A ODF Conn (m)']
f_conn_odf_b = interpolate.interp1d(tmp_x[~np.isnan(tmp_x)], tmp_y[~np.isnan(tmp_y)])
tmp_x = df['Conn R-Rsep omp (cm)']
tmp_y = df['A ODF Conn (m)']
f_conn_odf_c = interpolate.interp1d(tmp_x[~np.isnan(tmp_x)], tmp_y[~np.isnan(tmp_y)])

# The connections lengths at the same R-Rsep omp as the LP data points (and thus
# the sound speed data points).
cutoff_idx = 90
interp_conn_odf_a = f_conn_odf_a(df['A R-Rsep omp (cm)'][:-cutoff_idx])
interp_conn_odf_b = f_conn_odf_b(df['B R-Rsep omp (cm)'][:-cutoff_idx])
interp_conn_odf_c = f_conn_odf_c(df['C R-Rsep omp (cm)'][:-cutoff_idx])

# Now do the same as above just for idf side.
tmp_x = df['Conn R-Rsep omp (cm)']
tmp_y = df['A IDF Conn (m)']
f_conn_idf_a = interpolate.interp1d(tmp_x[~np.isnan(tmp_x)], tmp_y[~np.isnan(tmp_y)])
tmp_x = df['Conn R-Rsep omp (cm)']
tmp_y = df['A IDF Conn (m)']
f_conn_idf_b = interpolate.interp1d(tmp_x[~np.isnan(tmp_x)], tmp_y[~np.isnan(tmp_y)])
tmp_x = df['Conn R-Rsep omp (cm)']
tmp_y = df['A IDF Conn (m)']
f_conn_idf_c = interpolate.interp1d(tmp_x[~np.isnan(tmp_x)], tmp_y[~np.isnan(tmp_y)])
interp_conn_idf_a = f_conn_odf_a(df['A R-Rsep omp (cm)'][:-cutoff_idx])
interp_conn_idf_b = f_conn_odf_b(df['B R-Rsep omp (cm)'][:-cutoff_idx])
interp_conn_idf_c = f_conn_odf_c(df['C R-Rsep omp (cm)'][:-cutoff_idx])

# Dataframe to hold the answers we want.
dperp_df = pd.DataFrame()
dperp_df['A R-Rsep omp (cm)'] = df['A R-Rsep omp (cm)'][:-cutoff_idx]
dperp_df['B R-Rsep omp (cm)'] = df['B R-Rsep omp (cm)'][:-cutoff_idx]
dperp_df['C R-Rsep omp (cm)'] = df['C R-Rsep omp (cm)'][:-cutoff_idx]
dperp_df['A ODF Conn (m)']    = interp_conn_odf_a
dperp_df['B ODF Conn (m)']    = interp_conn_odf_b
dperp_df['C ODF Conn (m)']    = interp_conn_odf_c
dperp_df['A IDF Conn (m)']    = interp_conn_idf_a
dperp_df['B IDF Conn (m)']    = interp_conn_idf_b
dperp_df['C IDF Conn (m)']    = interp_conn_idf_c
dperp_df['Sound Speed (m/s)'] = df['Sound Speed (m/s)'][:-cutoff_idx]

# The six lambdas in omp coord. (m): [AODF, AIDF, BODF, BIDF, CODF, CIDF]
lambdas = np.array([0.0175, 0.0083, 0.023, 0.01, 0.013, 0.0063])

# Calculate Dperp.
dperp_df['A ODF Dperp (m2/s)'] = lambdas[0]**2 * dperp_df['Sound Speed (m/s)'] / dperp_df['A ODF Conn (m)']
dperp_df['A IDF Dperp (m2/s)'] = lambdas[1]**2 * dperp_df['Sound Speed (m/s)'] / dperp_df['A IDF Conn (m)']
dperp_df['B ODF Dperp (m2/s)'] = lambdas[2]**2 * dperp_df['Sound Speed (m/s)'] / dperp_df['B ODF Conn (m)']
dperp_df['B IDF Dperp (m2/s)'] = lambdas[3]**2 * dperp_df['Sound Speed (m/s)'] / dperp_df['B IDF Conn (m)']
dperp_df['C ODF Dperp (m2/s)'] = lambdas[4]**2 * dperp_df['Sound Speed (m/s)'] / dperp_df['C ODF Conn (m)']
dperp_df['C IDF Dperp (m2/s)'] = lambdas[5]**2 * dperp_df['Sound Speed (m/s)'] / dperp_df['C IDF Conn (m)']

# A normal plot.
if True:
    # Plot ODF and IDF on same graph.
    plt.rcParams.update({'font.size': 36})
    plt.plot(dperp_df['A R-Rsep omp (cm)'], dperp_df['A ODF Dperp (m2/s)'], 'C1-',  label='A-OTF')
    plt.plot(dperp_df['B R-Rsep omp (cm)'], dperp_df['B ODF Dperp (m2/s)'], 'C2-',  label='B-OTF')
    plt.plot(dperp_df['C R-Rsep omp (cm)'], dperp_df['C ODF Dperp (m2/s)'], 'C3-',  label='C-OTF')
    plt.plot(dperp_df['A R-Rsep omp (cm)'], dperp_df['A IDF Dperp (m2/s)'], 'C1--', label='A-ITF')
    plt.plot(dperp_df['B R-Rsep omp (cm)'], dperp_df['B IDF Dperp (m2/s)'], 'C2--', label='B-ITF')
    plt.plot(dperp_df['C R-Rsep omp (cm)'], dperp_df['C IDF Dperp (m2/s)'], 'C3--', label='C-ITF')
    plt.xlabel('R-Rsep omp (cm)')
    plt.ylabel(r'$\mathrm{D_{\perp}\ (m^2/s)}$')
    plt.legend(prop={"size":26})
    plt.show()

# A semilogy plot.
if True:
    plt.rcParams.update({'font.size': 36})
    plt.semilogy(dperp_df['A R-Rsep omp (cm)'], dperp_df['A ODF Dperp (m2/s)'], 'C1-',  label='A-OTF')
    plt.semilogy(dperp_df['B R-Rsep omp (cm)'], dperp_df['B ODF Dperp (m2/s)'], 'C2-',  label='B-OTF')
    plt.semilogy(dperp_df['C R-Rsep omp (cm)'], dperp_df['C ODF Dperp (m2/s)'], 'C3-',  label='C-OTF')
    plt.semilogy(dperp_df['A R-Rsep omp (cm)'], dperp_df['A IDF Dperp (m2/s)'], 'C1--', label='A-ITF')
    plt.semilogy(dperp_df['B R-Rsep omp (cm)'], dperp_df['B IDF Dperp (m2/s)'], 'C2--', label='B-ITF')
    plt.semilogy(dperp_df['C R-Rsep omp (cm)'], dperp_df['C IDF Dperp (m2/s)'], 'C3--', label='C-ITF')
    plt.xlabel('R-Rsep omp (cm)')
    plt.ylabel(r'$\mathrm{D_{\perp}\ (m^2/s)}$')
    plt.legend(prop={"size":26})
    plt.show()
