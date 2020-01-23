import get_lp
import get_ts
from scipy.signal import lfilter


# Grab ts_data.
ts_dict = get_ts.run_script(167196, 'core')

ts_psins = ts_dict['psins']['avg_psins']

# Grab lp data.
lp_dict = get_lp.get_dict_of_lps(167196)

for p in lp_dict.keys():
    min_idx = np.where(lp_dict[p]['time'] > 2500)[0][0]
    max_idx = np.where(lp_dict[p]['time'] > 5000)[0][0]
    avg_psin = lp_dict[p]['psin'][min_idx:max_idx].mean()
    lp_dict[p]['avg_psin'] = avg_psin

# probe 23 1.019 --> chord 8 1.022
# probe 25 1.028 --> chord 7 1.029
# probe 29 1.050 --> chord 5 1.056
# probe 31 1.063 --> chord 4 1.070
# probe 33 1.077 --> chord 4 1.070
# probe 35 1.092 --> chord 3 1.087

# Plot the time series of matching probes.
lp_min_idx = np.where(lp_dict['probe 25']['time'] > 2500)[0][0]
lp_max_idx = np.where(lp_dict['probe 25']['time'] > 5000)[0][0]
ts_min_idx = np.where(ts_dict['density']['X'] > 2500)[0][0]
ts_max_idx = np.where(ts_dict['density']['X'] > 5000)[0][0]
lp25_time = lp_dict['probe 25']['time'][lp_min_idx:lp_max_idx]
lp25_jsat = lp_dict['probe 25']['jsat'][lp_min_idx:lp_max_idx]
ts07_time = ts_dict['density']['X'][ts_min_idx:ts_max_idx]
ts07_dens = ts_dict['density']['Y'][6][ts_min_idx:ts_max_idx] * 10**(-18)

# Get the relative fluctation levels to the std dev. LP looks to sample at around 1 kHz. TS around 0.2 kHz.
lp25_jsat = (lp25_jsat - lp25_jsat.mean()) / lp25_jsat.std()
ts07_dens = (ts07_dens - ts07_dens.mean()) / ts07_dens.std()

# Filter jsat doing whatever this does.
#n = 10
#b = [1.0 / n] * n
#a = 1
#lp25_jsat_filt = lfilter(b, a, lp25_jsat)
fc = 0.4
fs = 1
w = fc / (fs / 2)
b, a = signal.butter(5, w, 'low')
lp25_jsat_filt = signal.filtfilt(b, a, lp25_jsat)

fig, ax = plt.subplots()
ax.plot(lp25_time, lp25_jsat_filt, label='Probe 25')
ax.plot(ts07_time, ts07_dens, label='Chord 7')
