import oedge_plots
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# Load the RCP data. Data has already been shifted inward by 1.5 cm due to EFIT uncertainties. Removing a couple
# bad data points.
rcp_path = "/Users/zamperini/Documents/d3d_work/divimp_files/oedge_tutorial/rcp_156195_2.csv"
rcp = pd.read_csv(rcp_path).iloc[:-4]

# Load OEDGE runs from F_fric scan, pull profile of ne at the RCP location, store in dictionary.
op_root = "/Users/zamperini/Documents/d3d_work/divimp_files/oedge_tutorial/"
ne_profs = {}
frics = np.arange(0.05, 1.00, 0.05)
for i in range(1, 20):
    op_path = "{}d3d-167196-osm-v1-mom{}.nc".format(op_root, i)
    op = oedge_plots.OedgePlots(op_path)
    ne_profs[frics[i-1]] = op.along_line(2.18, 2.30, -0.188, -0.188, "KNBS", "psin")

# For each ring, create an interpolation function of F_fric vs ne@RCP if possible.
f_f = {}
for ir in range(0, op.nrs):
    ne_at_rcp = []
    for fric in frics:

        # Mask for this ring. ir+1 is because OEDGE rings are 1-indexed, python is 0-indexed
        mask = np.array(ne_profs[fric]["ring"]) == ir+1

        # This shouldn't happen, but I (Shawn) haven't figured it out yet. It seems to not matter too much.
        if mask.sum() > 1:
            print("Warning! More than one value for ring {}".format(ir+1))
        if mask.sum() == 1:
            ne_at_rcp.append(float(np.array(ne_profs[fric]["KNBS"])[mask]))

    if len(ne_at_rcp) == 0:
        continue
    else:

        # Create interpolation function for F_fric(ne) so we can see what F_fric is needed for a desired ne value at
        # the location of the RCP.
        f_f[ir+1] = interp1d(ne_at_rcp, frics)

# For each ring with an interpolation function of F_fric(ne@RCP) find out what F_fric is needed to reproduce
# the RCP measurements. To do this we need an interpolation function of RCP_ne(psin).
f_rcp_ne = interp1d(rcp["psin"], rcp["Ne(E18 m-3)"] * 1e18)
fric_needed = {}
for ir in f_f.keys():

    # Get the ring's psin value so we can plug it into f_rcp_ne and get the desired density from OEDGE at the
    # RCP location.
    ring_psin = op.nc["PSIFL"][ir-1][0]  # 1-indexed to 0-indexed
    try:
        rcp_ne = f_rcp_ne(ring_psin)
        fric_needed[ir] = f_f[ir](rcp_ne)
    except ValueError:
        print("Ring {}: Outside of RCP data range - no value for F_fric given".format(ir))

# Now print out the data in a format that can be copy/pasted into input option *282.
print("'+242  Friction factor for momentum loss formula          '  1.0        # Default behavior is no momentum losses")
print("'*282  Momentum loss - ring specification                 '")
print("' ' '  Momentum loss - ring specification (dummy line)    '")
print("'  Ring   Ffric1     L1  Ffric2      L2    Number of rows:'  {}          # Only these rings have momentum losses".format(len(fric_needed)))
for ring, fric in fric_needed.items():
    print("     {}     {:.2f}    0.1    {:.2f}     0.1                   ".format(ring, fric, fric))