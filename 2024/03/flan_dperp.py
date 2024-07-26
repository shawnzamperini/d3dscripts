import flan_plots
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, BSpline
import matplotlib.ticker as mticker


def calc_dperp(flan_path):
    fp = flan_plots.FlanPlots(flan_path)
    fp_data = fp.plot_profiles(["vx", "nz"], plot_z=0.3125, 
        normtype=["log","log"], vmin=[1e18, 1], vmax=[2e19, 1e3], 
        skip_video=True)
    x = fp_data["x"]

    # Impurity density and x velocity averaged over the y coordinate.
    vx = np.nanmean(fp_data["data"][0], axis=2)
    nz = fp_data["data"][1].mean(axis=2)

    # Pick a range of frames after some time to equlibriate.
    # fstart = 75  # This is what I've normally been using
    fstart = 300  # Potentially a more steady-state onwards time
    vx_avg = np.nanmean(vx[fstart:], axis=0)

    # Should probably not include 0's? Idk, leaving alone for now. Results
    # aren't overall too different besids a bit higher Dperp. 
    nz_avg = nz[fstart:].mean(axis=0)
    #nz[nz == 0] = np.nan
    #nz_avg = np.nanmean(nz[fstart:], axis=0)

    def spline_fit(x, y, s, log=False, npoints=100):
        
        # Mask to remove any nans that sneak in.
        mask = ~np.isnan(y)
        x = x[mask]
        y = y[mask]
        
        if log:
            tck = splrep(x, np.log(y), s=s)
            xnew = np.linspace(x.min(), x.max(), npoints)
            return xnew, np.exp(BSpline(*tck)(xnew))
        else:
            tck = splrep(x, y, s=s)
            xnew = np.linspace(x.min(), x.max(), npoints)
            return xnew, BSpline(*tck)(xnew)
            
    # Spline fit to smooth the data some.
    x_spl, nz_avg_spl = spline_fit(x, nz_avg, 50, log=True) 
    x_spl, vx_avg_spl = spline_fit(x, vx_avg, 5e9, log=False) 

    # Calculate the average radial diffusion coefficient.
    gamma_x = nz_avg_spl * vx_avg_spl
    dnzdx = np.gradient(nz_avg_spl, x_spl)
    diff_x = -gamma_x / dnzdx

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, figsize=(9, 4))

    ax1.plot(x, vx_avg)
    ax1.plot(x_spl, vx_avg_spl)
    #ax1.set_yscale("log")

    ax2.plot(x, nz_avg)
    ax2.plot(x_spl, nz_avg_spl)
    ax2.set_yscale("log")

    ax3.plot(x_spl, diff_x)
    ax3.set_yscale("log")
    ax3.set_ylim([1, 100])

    fig.tight_layout()
    fig.show()


    return x_spl, diff_x, nz_avg_spl, vx_avg_spl


#flan_path = "/Users/zamperini/flandir/reg_testcase1/saved_results_5/reg_testcase1.nc"
coll_x, coll_diff, coll_nz, coll_vx = calc_dperp("/Users/zamperini/flandir/coll_on/coll_on.nc")
nocoll_x, nocoll_diff, nocoll_nz, nocoll_vx = calc_dperp("/Users/zamperini/flandir/coll_off/coll_off.nc")

# A standalone plot of the radial diffusion coefficient for pretty purposes.
xoffset = -2.259
fontsize = 14
fig, ax1 = plt.subplots(figsize=(5, 4))

# Shaded region of starting location.
xstart = 2.3025 + xoffset
#ax1.axvspan(xstart-0.001, xstart+0.001, color="tab:red", alpha=0.3)

ax1.plot(coll_x + xoffset, coll_diff, lw=4, color="k")
ax1.plot(coll_x + xoffset, coll_diff, lw=3, color="tab:red", label="ON")

# Comment out these three lines to just show the radial diffusion coefficient
# of the collision simulation.
#ax1.plot(nocoll_x + xoffset, nocoll_diff, lw=4, color="k")
#ax1.plot(nocoll_x + xoffset, nocoll_diff, lw=3, color="tab:purple", label="OFF")
#ax1.legend(fontsize=fontsize)

ax1.set_yscale("log")
ax1.set_ylim([1, 300])
ax1.set_xlim([2.31 + xoffset, None])
ax1.grid(which="both", alpha=0.2)
ax1.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{D_{r}\ (m^2/s)}$", fontsize=fontsize)
ax1.tick_params(axis='both', which='major', labelsize=fontsize - 2)
ax1.yaxis.set_major_formatter(mticker.ScalarFormatter())
fig.tight_layout()
fig.show()

# Plot of the density.
fig, ax1 = plt.subplots(figsize=(5, 4))

# Shaded region of starting location.
xstart = 2.3025 + xoffset
#ax1.axvspan(xstart-0.001, xstart+0.001, color="tab:red", alpha=0.3)

ax1.plot(coll_x + xoffset, coll_nz, lw=4, color="k")
ax1.plot(coll_x + xoffset, coll_nz, lw=3, color="tab:red", label="ON")

# Comment out these three lines to just show the radial diffusion coefficient
# of the collision simulation.
ax1.plot(nocoll_x + xoffset, nocoll_nz, lw=4, color="k")
ax1.plot(nocoll_x + xoffset, nocoll_nz, lw=3, color="tab:purple", label="OFF")
ax1.legend(fontsize=fontsize)

ax1.set_yscale("log")
#ax1.set_ylim([1, 300])
ax1.set_xlim([2.31 + xoffset, None])
ax1.grid(which="both", alpha=0.2)
ax1.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{n_{W}\ (a.u.)}$", fontsize=fontsize)
ax1.tick_params(axis='both', which='major', labelsize=fontsize - 2)
fig.tight_layout()
fig.show()


# Plot of the radial velocity.
fig, ax1 = plt.subplots(figsize=(5, 4))

# Shaded region of starting location.
xstart = 2.3025 + xoffset
#ax1.axvspan(xstart-0.001, xstart+0.001, color="tab:red", alpha=0.3)

ax1.plot(coll_x + xoffset, coll_vx, lw=4, color="k")
ax1.plot(coll_x + xoffset, coll_vx, lw=3, color="tab:red", label="ON")

#ax1.set_yscale("log")
#ax1.set_ylim([1, 300])
#ax1.set_xlim([2.31 + xoffset, None])
ax1.grid(which="both", alpha=0.2)
ax1.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{v_{r}\ (m/s)}$", fontsize=fontsize)
ax1.tick_params(axis='both', which='major', labelsize=fontsize - 2)
fig.tight_layout()
fig.show()



