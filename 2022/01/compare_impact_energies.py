import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerBase


casea_path = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-019a.nc"
casee_path = "/Users/zamperini/Documents/d3d_work/sput_testing/167196-sput-022a.nc"
casea = netCDF4.Dataset(casea_path)
casee = netCDF4.Dataset(casee_path)


class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle,
                       x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0,y0+width], [0.7*height,0.7*height],
                           linestyle=orig_handle[1], color='tab:red')
        l2 = plt.Line2D([x0,y0+width], [0.3*height,0.3*height],
                           linestyle=orig_handle[1], color='tab:purple')
        return [l1, l2]


def get_impact_flux(nc):

    # Load relevant variables.
    odouts = nc.variables["ODOUTS"][:]
    qtembs = nc.variables["QTEMBS"][:][0][:-2]  # Last 2 values seem to be junk.
    qtembsi = nc.variables["QTEMBSI"][:][0][:-2]
    qrnbs = nc.variables["QRNBS"][:][0][:-2]

    # The data is only for a single side, so mirror it.
    half = int(len(qtembs)/2)
    qtembs = np.append(qtembs[half:], qtembs[half:][::-1])
    qtembsi = np.append(qtembsi[half:], qtembsi[half:][::-1])
    qrnbs = np.append(qrnbs[half:], qrnbs[half:][::-1])


    # Calculate Eimpact and D flux, toroidal angle included.
    eimp = 2 * qtembsi + 3 * qtembs
    cs = 9.79e3 * np.sqrt((qtembs + qtembsi) / 2)
    csintb = np.sin(np.radians(14))
    flux = qrnbs * cs * csintb

    return odouts, eimp, flux


xa, eimpa, fluxa = get_impact_flux(casea)
xb, eimpb, fluxb = get_impact_flux(casee)

# Set equal number of ticks so grid lines up.
eimp_ticks = [0, 30, 60, 90, 120]
flux_ticks = [0, 3e22, 6e22, 9e22, 12e22]

fig, ax1 = plt.subplots(figsize=(5,4))

ax11 = ax1.twinx()

ax1.plot(xa, eimpa, color="tab:red")
ax1.plot(xb, eimpb, color="tab:red", linestyle="--")
ax11.plot(xa, fluxa, color="tab:purple")
ax11.plot(xb, fluxb, color="tab:purple", linestyle="--")

ax1.set_yticks(eimp_ticks)
ax11.set_yticks(flux_ticks)
ax1.grid()
ax1.set_ylim([0, 150])
ax11.set_ylim(0, 1.5e23)

ax1.set_xlabel("Distance along limiter (m)")
ax1.set_ylabel("Impact Energy (eV)", color="tab:red")
ax11.set_ylabel("D Flux (m-2 s-1)", color="tab:purple")
ax1.tick_params(axis='y', labelcolor="tab:red")
ax11.tick_params(axis='y', labelcolor="tab:purple")

plt.legend([("r","-"), ("limegreen","--")], ['Best Case', "Worst Case"],
           handler_map={tuple: AnyObjectHandler()})

fig.tight_layout()
fig.show()
