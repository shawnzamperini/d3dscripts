import pretty_plots as pp
from gadata import gadata
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'DejaVu Sans'

def get_data(tag):
    """
    Get the data for the tag fo rboth shots, storing in a list.
    """
    ga_247 = gadata(tag, 167247)
    ga_277 = gadata(tag, 167277)
    return ga_247.xdata, ga_247.zdata, ga_277.xdata, ga_277.zdata

def plot_data(ax, data, ylabel, ymult=1.0):
    ax.plot(data[0], data[1]*ymult, color=pp.tableau20[8], lw=2)
    ax.plot(data[2], data[3]*ymult, color=pp.tableau20[6], lw=2)
    ax.set_ylabel(ylabel)
    ax.set_xlim([0, 6000])

dens = get_data('DENSV2')
rvsout = get_data('RVSOUT')
pinj = get_data('PINJ')
ip = get_data('IP')
neped = get_data('PRMTAN_NEPED')
teped = get_data('PRMTAN_TEPED')
fs04 = get_data('FS04')

fig = plt.figure(figsize=(10,8))

ax1 = fig.add_subplot(421)
plot_data(ax1, dens, 'Density (m-3)')

ax2 = fig.add_subplot(425, sharex=ax1)
plot_data(ax2, rvsout, 'Strike Point (m)')

ax3 = fig.add_subplot(423, sharex=ax1)
plot_data(ax3, pinj, 'Pinj (MW)', ymult=0.001)

ax4 = fig.add_subplot(427, sharex=ax1)
plot_data(ax4, ip, 'IP (MA)', ymult=1e-6)

ax5 = fig.add_subplot(422, sharex=ax1)
plot_data(ax5, neped, 'Pedestal ne (m-3)')

ax6 = fig.add_subplot(424, sharex=ax1)
plot_data(ax6, teped, 'Pedestal Te (eV)')
ax6.set_xlabel('Time (ms)')

ax7 = fig.add_subplot(426, sharex=ax1)
plot_data(ax7, fs04, 'FS04 (ph/cm2/sr/s)')

fig.tight_layout()
fig.show()
