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
f_conn_idf_a = interpolate.interp1d(df['conn_omp_a'], df['conn_idf_a'])
f_conn_idf_b = interpolate.interp1d(df['conn_omp_b'], df['conn_idf_b'])
f_conn_idf_c = interpolate.interp1d(df['conn_omp_c'], df['conn_idf_c'])
interp_conn_idf_a = f_conn_idf_a(df['samp_omp_a'])
interp_conn_idf_b = f_conn_idf_b(df['samp_omp_b'])
interp_conn_idf_c = f_conn_idf_c(df['samp_omp_c'])

coll_omp_a = df['samp_omp_a']
coll_omp_b = df['samp_omp_b']
coll_omp_c = df['samp_omp_c']
coll_odf_a = np.zeros(len(coll_omp_a))
coll_odf_b = np.zeros(len(coll_omp_b))
coll_odf_c = np.zeros(len(coll_omp_c))
coll_idf_a = np.zeros(len(coll_omp_a))
coll_idf_b = np.zeros(len(coll_omp_b))
coll_idf_c = np.zeros(len(coll_omp_c))

# Get A ODF sampling length (mislabelled as coll).
idx = 0
for samp, odf in zip(df['samp_l_a'], interp_conn_odf_a):
    if np.isnan(odf):
        coll_odf_a[idx] = np.nan
    else:
        tmp_coll = min([samp, odf])
        #print(str(samp) + ' : ' + str(conn) + ' --> ' + str(tmp_coll))
        coll_odf_a[idx] = tmp_coll
    idx += 1
# Get B ODF sampling length.
idx = 0
for samp, odf in zip(df['samp_l_b'], interp_conn_odf_b):
    if np.isnan(odf):
        coll_odf_b[idx] = np.nan
    else:
        tmp_coll = min([samp, odf])
        #print(str(samp) + ' : ' + str(conn) + ' --> ' + str(tmp_coll))
        coll_odf_b[idx] = tmp_coll
    idx += 1
# Get C ODF sampling length.
idx = 0
for samp, odf, in zip(df['samp_l_c'], interp_conn_odf_c):
    if np.isnan(odf):
        coll_odf_c[idx] = np.nan
    else:
        tmp_coll = min([samp, odf])
        #print(str(samp) + ' : ' + str(conn) + ' --> ' + str(tmp_coll))
        coll_odf_c[idx] = tmp_coll
    idx += 1



# Get A IDF sampling length.
idx = 0
for samp, idf in zip(df['samp_l_a'], interp_conn_idf_a):
    if np.isnan(idf):
        coll_idf_a[idx] = np.nan
    else:
        tmp_coll = min([samp, idf])
        #print(str(samp) + ' : ' + str(conn) + ' --> ' + str(tmp_coll))
        coll_idf_a[idx] = tmp_coll
    idx += 1
# Get A IDF sampling length.
idx = 0
for samp, idf in zip(df['samp_l_b'], interp_conn_idf_b):
    if np.isnan(idf):
        coll_idf_b[idx] = np.nan
    else:
        tmp_coll = min([samp, idf])
        #print(str(samp) + ' : ' + str(conn) + ' --> ' + str(tmp_coll))
        coll_idf_b[idx] = tmp_coll
    idx += 1
# Get A IDF sampling length.
idx = 0
for samp, idf in zip(df['samp_l_c'], interp_conn_idf_c):
    if np.isnan(idf):
        coll_idf_c[idx] = np.nan
    else:
        tmp_coll = min([samp, idf])
        #print(str(samp) + ' : ' + str(conn) + ' --> ' + str(tmp_coll))
        coll_idf_c[idx] = tmp_coll
    idx += 1

linewidth = 5
if True:
    # Graph 1. Note coll is actually the sampling length and vice versa.
    plt.rcParams.update({'font.size': 42})
    plt.rcParams.update({'figure.autolayout': True})
    plt.semilogy(df['samp_omp_a'], df['samp_l_a'],   '-', linewidth=linewidth, label='A Collection length')
    plt.semilogy(df['conn_omp_a'], df['conn_odf_a'], '--', linewidth=linewidth, label='OTF Connection length')
    plt.semilogy(df['conn_omp_a'], df['conn_idf_a'], '--', linewidth=linewidth, label='ITF Connection length')
    #plt.semilogy(coll_omp_a,       coll_a,              label='Sampling length')
    plt.legend(prop={"size":26})
    plt.xlim([0, 17.5])
    plt.xlabel("R-Rsep omp (cm)")
    plt.ylabel("Length (m)")
    plt.title("Probe Lengths")
    #plt.savefig('coll_otf_itf.png', bbox_inches='tight')
    plt.show()

if True:
    # Graph 2
    plt.rcParams.update({'font.size': 42})
    plt.semilogy(coll_omp_a, coll_odf_a, linewidth=linewidth, label='A')
    plt.semilogy(coll_omp_b, coll_odf_b, linewidth=linewidth, label='B')
    plt.semilogy(coll_omp_c, coll_odf_c, linewidth=linewidth, label='C')
    plt.legend(prop={"size":26})
    plt.xlim([0, 17.5])
    plt.xlabel("R-Rsep omp (cm)")
    plt.ylabel("Length (m)")
    plt.title("OTF Sampling Lengths")
    plt.savefig('samp_otf.png')
    plt.show()

if True:
    # Graph 3. Same as above except IDF.
    plt.rcParams.update({'font.size': 42})
    plt.semilogy(coll_omp_a, coll_idf_a, linewidth=linewidth, label='A')
    plt.semilogy(coll_omp_b, coll_idf_b, linewidth=linewidth, label='B')
    plt.semilogy(coll_omp_c, coll_idf_c, linewidth=linewidth, label='C')
    plt.legend(prop={"size":26})
    plt.xlim([0, 17.5])
    plt.xlabel("R-Rsep omp (cm)")
    plt.ylabel("Length (m)")
    plt.title("ITF Sampling Lengths")
    plt.savefig('samp_itf.png')
    plt.show()

if True:
    plt.rcParams.update({'font.size': 42})
    plt.semilogy(coll_omp_a, coll_odf_a, color='C0', linewidth=linewidth, label='A OTF')
    plt.semilogy(coll_omp_b, coll_odf_b, color='C1', linewidth=linewidth, label='B OTF')
    plt.semilogy(coll_omp_c, coll_odf_c, color='C2', linewidth=linewidth, label='C OTF')
    plt.semilogy(coll_omp_a, coll_idf_a, '--', color='C0', linewidth=linewidth, label='A ITF')
    plt.semilogy(coll_omp_b, coll_idf_b, '--', color='C1', linewidth=linewidth, label='B ITF')
    plt.semilogy(coll_omp_c, coll_idf_c, '--', color='C2', linewidth=linewidth, label='C ITF')
    plt.legend(prop={"size":26})
    plt.xlim([0, 17.5])
    plt.xlabel("R-Rsep omp (cm)")
    plt.ylabel("Length (m)")
    plt.savefig('samp_otf_itf.png')
    plt.show()
