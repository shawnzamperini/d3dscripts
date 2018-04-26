# IPython log file

get_ipython().magic('run ts_comparision.py')
get_ipython().magic('quickref')
get_ipython().magic('logstart')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
get_ipython().magic('matplotlib')
def myplot():
    ms=7
    cs=5
    alpha=0.5
    fig = plt.figure()
    ax1 = plt.subplot(111)
    params = {
       'axes.labelsize': 8,
       'font.size': 8,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha)
    params = {
       'axes.labelsize': 8,
       'font.size': 8,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    params = {
       'axes.labelsize': 8,
       'font.size': 8,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    params = {
       'axes.labelsize': 8,
       'font.size': 8,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    ax1.set_xlim(0.95, 1.35)
    ax1.set_ylim(0, 80)
    params = {
       'axes.labelsize': 8,
       'font.size': 8,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    ax1.set_xlim(1, 1.35)
    ax1.set_ylim(0, 80)
    params = {
       'axes.labelsize': 8,
       'font.size': 8,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$\mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    ax1.set_xlim(1, 1.35)
    ax1.set_ylim(0, 80)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$\mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    ax1.set_xlim(1, 1.35)
    ax1.set_ylim(0, 80)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$\mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    ax1.set_xlim(1, 1.35)
    ax1.set_ylim(0, 80)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
      # 'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
myplot()
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$\mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    ax1.set_xlim(1, 1.35)
    ax1.set_ylim(0, 80)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$\mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    ax1.set_xlim(1, 1.35)
    ax1.set_ylim(0, 80)
    params = {
       'axes.labelsize': 14,
       'font.size': 14,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$\mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    ax1.set_xlim(1, 1.35)
    ax1.set_ylim(0, 80)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
myplot()
def myplot():
    ms=7
    cs=5
    alpha=0.5
    st=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.errorbar(x=psin_lp_195p2[::st], y=Te_lp_195p2[::st],
                 xerr=0.01, yerr=(Te_lp_195p2*0.2)[::st],
                 capsize=cs, ms=ms, fmt='r^', alpha=alpha, label='Langmuir Probe')
    ax2 = plt.subplot(111)
    ax2.errorbar(x=ts196['psins']['avg_psins'], y=ts192['psins']['avg_Tes'],
                 xerr=0.01, yerr=ts192['psins']['avg_Tes_err'], capsize=cs,
                 ms=ms, fmt='b^', alpha=alpha, label='Thomson Scattering')
    ax1.set_xlabel(r'$\mathrm{\phi_n}$')
    ax1.set_ylabel('Te (eV)')
    ax1.legend()
    ax1.set_xlim(1, 1.35)
    ax1.set_ylim(0, 80)
    params = {
       'axes.labelsize': 16,
       'font.size': 16,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
myplot()
myplot()
# Now creating plots of sampling lengths.
get_ipython().magic('run samp_coll_conn.py')
get_ipython().magic('run samp_coll_conn.py')
get_ipython().magic('run samp_coll_conn.py')
def ourplot():
    fig = plt.figure()
    ax1 = plt.subfigure(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, color='r')
    ax1.semilogy(coll_omp_b, coll_odf_b, color='b')
    ax1.semilogy(coll_omp_c, coll_odf_c, color='g')
    ax1.semilogy(coll_omp_a, coll_idf_a, color='r')
    ax1.semilogy(coll_omp_b, coll_idf_b, color='b')
    ax1.semilogy(coll_omp_c, coll_idf_c, color='g')
    
ourplot()
coll_omp_a
coll_odf_a
ourplot()
def ourplot():
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, color='r')
    ax1.semilogy(coll_omp_b, coll_odf_b, color='b')
    ax1.semilogy(coll_omp_c, coll_odf_c, color='g')
    ax1.semilogy(coll_omp_a, coll_idf_a, color='r')
    ax1.semilogy(coll_omp_b, coll_idf_b, color='b')
    ax1.semilogy(coll_omp_c, coll_idf_c, color='g')
    
outplot()
ourplot()
def ourplot():
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, fmt='r')
    ax1.semilogy(coll_omp_b, coll_odf_b, fmt='b')
    ax1.semilogy(coll_omp_c, coll_odf_c, fmt='g')
    ax1.semilogy(coll_omp_a, coll_idf_a, fmt='r')
    ax1.semilogy(coll_omp_b, coll_idf_b, fmt='b')
    ax1.semilogy(coll_omp_c, coll_idf_c, fmt='g')
    
ourplot()
def ourplot():
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, fmt='r-')
    ax1.semilogy(coll_omp_b, coll_odf_b, fmt='b-')
    ax1.semilogy(coll_omp_c, coll_odf_c, fmt='g-')
    ax1.semilogy(coll_omp_a, coll_idf_a, fmt='r--')
    ax1.semilogy(coll_omp_b, coll_idf_b, fmt='b--')
    ax1.semilogy(coll_omp_c, coll_idf_c, fmt='g--')
    
ourplot()
def ourplot():
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--')
    
ourplot()
def ourplot():
    lw=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--')
    params = {
       'axes.labelsize': 12,
       'text.fontsize': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
ourplot()
def ourplot():
    lw=2
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', linewidth=lw)
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', linewidth=lw)
    ax1.set_xlabel(r'$\mathrm{R - R_{sep} omp (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.5
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw)
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.5
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw)
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.5
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw)
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    fig.tight_layout()
    
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.5
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw)
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw)
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw)
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw)
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    fig.tight_layout()
    
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.5
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    fig.tight_layout()
    
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.5
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both')
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    fig.tight_layout()
    
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.5
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    fig.tight_layout()
    
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 8,
       'font.size': 8,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    fig.tight_layout()
    
ourplot()
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 10,
       'font.size': 10,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    fig.tight_layout()
    
outplot()
ourplot()
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 10,
       'font.size': 10,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    #fig.tight_layout()
    
ourplot()
ourplot()
ourplot()
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    #fig.tight_layout()
    
ourplot()
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 14,
       'font.size': 14,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    #fig.tight_layout()
    
ourplot()
ourplot()
def ourplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 16,
       'font.size': 16,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    #fig.tight_layout()
    
ourplot()
ourplot()
# Now do one comparing the coll, conn and samp of A
def yourplot():
    lw=2
    alpha=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='Collection Length')
    ax1.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax1.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
yourplot()
def yourplot():
    lw=2
    alpha=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax1.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax1.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
yourplot()
yourplot()
def yourplot():
    lw=2
    alpha=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax1.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax1.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 14,
       'font.size': 14,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
yourplot()
yourplot()
def yourplot():
    lw=2
    alpha=1.0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax1.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax1.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    params = {
       'axes.labelsize': 16,
       'font.size': 16,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
yourplot()
yourplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(112)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    plt.grid(True, which='both', alpha=0.5)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    ax2 = plt.subplot(211)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    plt.grid(True, which='both', alpha=0.5)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    plt.grid(True, which='both', alpha=0.5)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    plt.grid(True, which='both', alpha=0.5)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    #plt.grid(True, which='both', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    plt.grid(True, which='both', alpha=0.5)
    
weplot()
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    #plt.grid(True, which='both', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    #plt.grid(True, which='both', alpha=0.5)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    #plt.grid(True, which='both', alpha=0.5)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    plt.grid(True, which='both', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='both', alpha=0.5)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='both', alpha=0.5)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    ax1.grid(True, which='major', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='major', alpha=0.5)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    ax1.grid(True, which='major', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='major', alpha=0.5)
    fig.tight_layout()
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    ax1.grid(True, which='major', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='major', alpha=0.5)
    fig.tight_layout()
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    ax1.grid(True, which='major', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='major', alpha=0.5)
    fig.tight_layout()
    params = {
       'axes.labelsize': 16,
       'font.size': 16,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
weplot()
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    ax1.grid(True, which='major', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    #ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='major', alpha=0.5)
    fig.tight_layout()
    params = {
       'axes.labelsize': 16,
       'font.size': 16,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    ax1.grid(True, which='major', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    #ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='major', alpha=0.5)
    fig.tight_layout()
    params = {
       'axes.labelsize': 14,
       'font.size': 14,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    ax1.grid(True, which='major', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    #ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='major', alpha=0.5)
    fig.tight_layout()
    params = {
       'axes.labelsize': 14,
       'font.size': 14,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       #'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
weplot()
weplot()
def weplot():
    lw=2
    alpha1 = 0.4
    alpha2=1.0
    fig, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
    #ax1 = plt.subplot(211)
    ax1.semilogy(coll_omp_a, coll_odf_a, 'r-', alpha=alpha1, linewidth=lw, label='A OTF')
    ax1.semilogy(coll_omp_b, coll_odf_b, 'b-', alpha=alpha1, linewidth=lw, label='B OTF')
    ax1.semilogy(coll_omp_c, coll_odf_c, 'g-', alpha=alpha1, linewidth=lw, label='C OTF')
    ax1.semilogy(coll_omp_a, coll_idf_a, 'r--', alpha=alpha2, linewidth=lw, label='A ITF')
    ax1.semilogy(coll_omp_b, coll_idf_b, 'b--', alpha=alpha2, linewidth=lw, label='B ITF')
    ax1.semilogy(coll_omp_c, coll_idf_c, 'g--', alpha=alpha2, linewidth=lw, label='C ITF')
    ax1.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax1.set_ylabel('Sampling Length (m)')
    ax1.legend()
    ax1.grid(True, which='major', alpha=0.5)
    #ax2 = plt.subplot(111)
    ax2.semilogy(df['samp_omp_a'], df['samp_l_a'], 'r-', linewidth=lw, 
                 label='A Collection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_odf_a'], 'g--', linewidth=lw, 
                 label='OTF Connection Length')
    ax2.semilogy(df['conn_omp_a'], df['conn_idf_a'], 'b--', linewidth=lw, 
                 label='ITF Connection Length')
    #ax2.set_xlabel(r'$\mathrm{R - R_{sep}\ omp\ (cm)}$')
    ax2.set_ylabel('Length (m)')
    ax2.legend()
    ax2.grid(True, which='major', alpha=0.5)
    fig.tight_layout()
    params = {
       'axes.labelsize': 12,
       'font.size': 12,
       'legend.fontsize': 10,
       'xtick.labelsize': 10,
       'ytick.labelsize': 10,
       'text.usetex': False,
       #'figure.figsize': [4.5, 4.5]
       }
    plt.rcParams.update(params)
    
weplot()
weplot()
weplot()
