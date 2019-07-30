from ThomsonClass import ThomsonClass
import numpy as np
import pretty_plots as pp


time_min = 4000
time_max = 5000

ts_div = ThomsonClass(167192, 'divertor')
ts_div.load_ts()

# Load in all the times.
div_times = ts_div.ts_dict['temp']['X'].round(2)

# Restrict to the time range input.
div_idx  = np.where(np.logical_and(div_times>time_min, div_times<time_max))
div_times  = div_times[div_idx]

#div_times = np.linspace(time_min,time_max,5)

# Map the data to EFIT.
ts_div.map_to_efit(times=div_times, trunc_div=True)

chord = 0
psin_bot = np.array([])
te_bot   = np.array([])
for time in div_times:
    coord = ts_div.temp_df_omp[str(time) + ' psin'][chord]
    psin_bot = np.append(psin_bot, coord[0])
    te_bot   = np.append(te_bot,   coord[1])

pp.pplot(psin_bot, te_bot, fmt='.')
