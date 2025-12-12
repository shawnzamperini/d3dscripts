from ThomsonClass import ThomsonClass
import numpy as np


# Inputs
# H-mode
#shot = 200083; tmin = 2000; tmax = 4000

# QH-mode
#shot = 201661; tmin = 2200; tmax = 3400

# Super H-mode
#shot = 203188; tmin = 2800; tmax = 3300

# High powered reverse L-mode
#shot = 200682; tmin = 2000; tmax = 5000

# From Lauren's experiment, just for fun
shot = 203795; tmin = 1500; tmax = 4000

# Negative triangularity
#shot = 193765; tmin = 5200; tmax = 5700

# High beta-p
#shot = 170361; tmin = 4000; tmax = 4600

tmany = 20

# EFIT04 is missing sometimes
if shot in [200682, 203795]:
	tree = "EFIT01"
else:
	tree = "EFIT04"

ts = ThomsonClass(shot, 'core')
ts.load_ts(tunnel=False)
ts.map_to_efit(np.linspace(tmin, tmax, tmany), tree=tree, trunc_div=False)

# Pull out the arrays.
rmrsomp  = ts.avg_omp['RminRsep_omp']
te = ts.avg_omp['Te_omp']
ne = ts.avg_omp['ne_omp']
rmrsomp_err  = ts.avg_omp['RminRsep_omp_err']
te_err = ts.avg_omp['Te_omp_err']
ne_err = ts.avg_omp['ne_omp_err']

