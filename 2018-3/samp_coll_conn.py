import pandas as pd
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

filename = '../Data/LPData.xlsx'
print('Loading Excel file\n  Filename: ' + filename)
df = pd.read_excel(filename, sheetname='Lengths',
                   usecols=[1,2,4,5,6,7,8,9,11,12,13,17,18,19,23,24,25],
                   skiprows=3,
                   names=['lp_r',       'lp_te',      'samp_omp_a', 'samp_omp_b',
                          'samp_omp_c', 'samp_l_a',   'samp_l_b',   'samp_l_c',
                          'conn_omp_a', 'conn_omp_b', 'conn_omp_c', 'conn_odf_a',
                          'conn_odf_b', 'conn_odf_c', 'conn_idf_a', 'conn_idf_b',
                          'conn_idf_c'])

f_conn_odf_a = interpolate.interp1d(df['conn_omp_a'], df['conn_odf_a'])
f_conn_odf_b = interpolate.interp1d(df['conn_omp_b'], df['conn_odf_b'])
f_conn_odf_c = interpolate.interp1d(df['conn_omp_c'], df['conn_odf_c'])
interp_conn_odf_a = f_conn_odf_a(df['samp_omp_a'])
interp_conn_odf_b = f_conn_odf_b(df['samp_omp_b'])
interp_conn_odf_c = f_conn_odf_c(df['samp_omp_c'])
f_conn_idf_a = interpolate.interp1d(df['conn_omp_a'], df['conn_odf_a'])
f_conn_idf_b = interpolate.interp1d(df['conn_omp_b'], df['conn_odf_b'])
f_conn_idf_c = interpolate.interp1d(df['conn_omp_c'], df['conn_odf_c'])
interp_conn_idf_a = f_conn_idf_a(df['samp_omp_a'])
interp_conn_idf_b = f_conn_idf_b(df['samp_omp_b'])
interp_conn_idf_c = f_conn_idf_c(df['samp_omp_c'])

coll_omp_a = df['samp_omp_a']
coll_omp_b = df['samp_omp_b']
coll_omp_c = df['samp_omp_c']
coll_a = np.zeros(len(coll_omp_a))
coll_b = np.zeros(len(coll_omp_b))
coll_c = np.zeros(len(coll_omp_c))

idx = 0
for samp, odf, idf in zip(df['samp_l_a'], interp_conn_odf_a, interp_conn_idf_a):
    if np.isnan(odf) or np.isnan(idf):
        coll_a[idx] = np.nan
    else:
        tmp_coll = min([samp, odf, idf])
        #print(str(samp) + ' : ' + str(conn) + ' --> ' + str(tmp_coll))
        coll_a[idx] = tmp_coll
    idx += 1

idx = 0
for samp, odf, idf in zip(df['samp_l_b'], interp_conn_odf_b, interp_conn_idf_b):
    if np.isnan(odf) or np.isnan(idf):
        coll_b[idx] = np.nan
    else:
        tmp_coll = min([samp, odf, idf])
        #print(str(samp) + ' : ' + str(conn) + ' --> ' + str(tmp_coll))
        coll_b[idx] = tmp_coll
    idx += 1

idx = 0
for samp, odf, idf in zip(df['samp_l_c'], interp_conn_odf_c, interp_conn_idf_c):
    if np.isnan(odf) or np.isnan(idf):
        coll_c[idx] = np.nan
    else:
        tmp_coll = min([samp, odf, idf])
        #print(str(samp) + ' : ' + str(conn) + ' --> ' + str(tmp_coll))
        coll_c[idx] = tmp_coll
    idx += 1

# Graph 1
plt.rcParams.update({'font.size': 34})
plt.semilogy(df['samp_omp_a'], df['samp_l_a'],   '--', label='Sampling length')
plt.semilogy(df['conn_omp_a'], df['conn_odf_a'], '--', label='ODF Connection length')
plt.semilogy(df['conn_omp_a'], df['conn_idf_a'], '--', label='IDF Connection length')
plt.semilogy(coll_omp_a,       coll_a,              label='Collection length')
plt.legend(prop={"size":22})
plt.xlabel("R-Rsep omp (cm)")
plt.ylabel("Length (m)")
plt.title("Probe Lengths")
plt.show()

# Graph 2
plt.rcParams.update({'font.size': 34})
plt.semilogy(coll_omp_a, coll_a, label='A')
plt.semilogy(coll_omp_b, coll_b, label='B')
plt.semilogy(coll_omp_c, coll_c, label='C')
plt.legend(prop={"size":22})
plt.xlabel("R-Rsep omp (cm)")
plt.ylabel("Length (m)")
plt.title("Collection Lengths")
plt.show()
