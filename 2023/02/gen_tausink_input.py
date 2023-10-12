# Script to help output some of the input for the 1DLIM model, I put the output into tausink_v2.xlsx.
import pandas as pd
from scipy.interpolate import interp1d


# Load in connection length data.
def load_mafot(itf_path, otf_path):
    print("Loading MAFOT data...")
    columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
               "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
    mafot_itf = pd.read_csv(itf_path, skiprows=52, names=columns, delimiter="\t")
    mafot_otf = pd.read_csv(otf_path, skiprows=52, names=columns, delimiter="\t")
    conns_r = mafot_itf["R (m)"]
    conns_l_itf = mafot_itf["Lconn (km)"].values * 1000  # km to m
    conns_l_otf = mafot_otf["Lconn (km)"].values * 1000
    return {"r":conns_r, "itf":conns_l_itf, "otf":conns_l_otf}


mafot = load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp2000_conn_-1.dat",
                   "/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp2000_conn_+1.dat")

# We already have R-Rsep coordinates setup that align with the R bins, so convert those to an R value for this shot.
rsep = 2.21602
df = pd.read_excel("/Users/zamperini/My Drive/Research/Documents/2023/02/tausink_v2.xlsx", sheet_name="167196_conns_2cm")
rs = df["r"].to_numpy() + rsep

# Setup interpolation function for MAFOT. Find the connection lengths at our coordinates.
fitf = interp1d(mafot["r"], mafot["itf"], bounds_error=False, fill_value=(mafot["itf"][0], mafot["itf"][-1]))
fotf = interp1d(mafot["r"], mafot["otf"], bounds_error=False, fill_value=(mafot["otf"][0], mafot["otf"][-1]))
conns = fitf(rs) + fotf(rs)
