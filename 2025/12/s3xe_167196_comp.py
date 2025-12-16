import numpy as np
import matplotlib.pyplot as plt
from S3X_EIRENE import S3X_EIRENE


# Load an S3XE case
simulation_folder = "/pscratch/sd/z/zamp/soledge3x_rundir/diiid_167196_v1"

sim = S3X_EIRENE()

sim.load_S3XE_simulation(
    base_folder             = simulation_folder,
    run_dir                 = "run_dir", # Default
    plasma_folder           = "Plasma", # Default
    load_restart            = False, # Default
    load_final              = False, # Default
    load_only_indices       = None, # Default, all plasma files. can be integer, or list/range (eg. range(0,1000,100)), put -1 to load only last file
    load_emergency_saves    = False, # Default
    force_reload            = False
)

# Load QuadPlasma object
qp = sim.SOLEDGE_plasma

# Load QuadFields
n_e = sim.SOLEDGE_plasma.species["e-"]["n"]
T_e = sim.SOLEDGE_plasma.species["e-"]["T"]
T_i = sim.SOLEDGE_plasma.species["D+"]["T"]

# Number of 
ne_RZ = qp.get_values_from_RZ(n_e, [[2.15, 0.0]], -1)
