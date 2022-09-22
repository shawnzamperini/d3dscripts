from gadata import gadata
import MDSplus as mds
import numpy as np
import sys


def main(params):

    # Pull out the info.
    shot = int(params[0])
    tmin = float(params[1])
    tmax = float(params[2])

    # Create MDSplus object and GA objs.
    conn = mds.Connection("atlas.gat.com")
    prad_obj = gadata("PRAD_CORE", shot, connection=conn)
    pinj_obj = gadata("PINJ", shot, connection=conn)
    pech_obj = gadata("ECHPWR", shot, connection=conn)
    pohm_obj = gadata("POH", shot, connection=conn)
    wmhd_obj = gadata("DWMHDF", shot, connection=conn)  # Change in stored energy

    # Calculate the average values.
    def calc_avg(obj, dbg_str):
        mask = np.logical_and(obj.xdata>tmin, obj.xdata<tmax)
        if len(obj.zdata[mask]) == 0:
            print("Warning: No data for {}".format(dbg_str))
            return 0
        else:
            return obj.zdata[mask].mean()
    avg_prad = calc_avg(prad_obj, "prad") / 1e6  # W to MW
    avg_pinj = calc_avg(pinj_obj, "pinj") / 1e3  # kW to MW
    avg_pech = calc_avg(pech_obj, "pech") / 1e6  # W to MW
    avg_pohm = calc_avg(pohm_obj, "pohm") / 1e6  # W to MW
    avg_wmhd = calc_avg(wmhd_obj, "wmhd") / 1e6  # W to MW

    # Finally return PSOL.
    psol = avg_pinj + avg_pech + avg_pohm + avg_wmhd - avg_prad
    print("PINJ = {:.4f}".format(avg_pinj))
    print("PRAD = {:.4f}".format(avg_prad))
    print("PECH = {:.4f}".format(avg_pech))
    print("POHM = {:.4f}".format(avg_pohm))
    print("PMHD = {:.4f}".format(avg_wmhd))
    print("PSOL = {:.4f}".format(psol))
    return psol


if __name__ == "__main__":
   main(sys.argv[1:])
