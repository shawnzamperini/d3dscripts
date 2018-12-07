from ThomsonClass import ThomsonClass
import numpy as np
import pretty_plots as pp
from scipy.optimize import curve_fit

def detach_front(shot, tmin=2000, tmax=4500, tmany=10, ref_time=2000,
                 detach_front_te=5, offset=0.01):
    ts = ThomsonClass(shot, 'divertor')
    ts.load_ts()
    ts.map_to_efit(times=np.linspace(tmin, tmax, tmany), ref_time=ref_time)
    ts.heatmap(detach_front_te=detach_front_te, offset=offset)

    return ts

def ts_falloff(shot, tmin, tmax, tmany):

    # Load the core TS.
    ts = ThomsonClass(shot, 'core')
    ts.load_ts()
    ts.map_to_efit(times=np.linspace(tmin, tmax, tmany), debug=False)

    try:
        avg_omps = np.array([])
        avg_tes  = np.array([])
        avg_nes  = np.array([])
        for chord in range(0, len(ts.temp_df_omp.index)):
            tmp_omps = np.array([])
            tmp_tes  = np.array([])
            tmp_nes  = np.array([])

            for time in range(0, tmany):

                # Get the tuple data point for this chord at this time (r-rsep_omp, Te).
                tmp_o = ts.temp_df_omp.values[chord][time][0]
                tmp_t = ts.temp_df_omp.values[chord][time][1]
                tmp_n = ts.dens_df_omp.values[chord][time][1]
                tmp_omps = np.append(tmp_omps, tmp_o)
                tmp_tes  = np.append(tmp_tes,  tmp_t)
                tmp_nes  = np.append(tmp_nes,  tmp_n)

            # Get the average for this chord. Append it to avg_omps/avg_tes.
            avg_omps = np.append(avg_omps, tmp_omps.mean())
            avg_tes  = np.append(avg_tes,  tmp_tes.mean())
            avg_nes  = np.append(avg_nes,  tmp_nes.mean())

        while True:
            # Plot it up and ask for where the pedestal starts. This is where the fit will start.
            print("Identify the start of the pedestal (q to quit).")
            fig = pp.pplot(avg_omps, avg_nes)
            try:
                ped_start = float(input("Start of pedestal: "))
            except:
                break

            # Find closest point to the entered value.
            ped_idx = np.where(np.abs(avg_omps-ped_start) == np.abs(avg_omps-ped_start).min())[0][0]

            # Start the fit from the first point to the right.
            x = avg_omps[:ped_idx+1]
            y = avg_nes[:ped_idx+1] * 1e-18

            def exp_fit(x, a, b):
                return a * np.exp(-x*b)

            popt_ne, pcov_ne = curve_fit(exp_fit, x, y-y.min(), maxfev=5000)
            x_fit = np.linspace(x.min(), x.max(), 100)
            y_fit = exp_fit(x_fit, *popt_ne) + y.min()

            fig = pp.pplot(x, y*1e18, xlabel='R-Rsep OMP (m)', ylabel='Density (m-3)')
            fig = pp.pplot(x_fit, y_fit*1e18, fmt='--', fig=fig, xlabel='R-Rsep OMP (m)',
                           ylabel='Density (m-3)')

            print("Density falloff length (cm): {:.2f}".format(100*1/popt_ne[1]))

        return ts

    except:
        return ts
