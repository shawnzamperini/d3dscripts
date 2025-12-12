import pickle
import pandas as pd
import numpy as np

# Load pickled gfile
"""
import pickle
import numpy as np
from os.path import expanduser
import omfit_classes

gfile = OMFIT["EFIT"]["FILES"]["gEQDSK"]
shot = int(gfile["CASE"][3].split("#")[1])
time = int(gfile["CASE"][4].split("ms")[0])
home = expanduser("~")
fname = "{}/{}_{}.pickle".format(home, shot, time)

with open(fname, "wb") as f:
    gfile_out = {}
    for k, v in gfile.items():
        print("{} {}".format(k, type(v)))
        if k in ["AuxNamelist", "fluxSurfaces", "_desc"]:
            continue
        if type(v) == omfit_classes.sortedDict.SortedDict:
            for k2, v2 in v.items():
                print("  {} {}".format(k2, type(v2)))
                gfile_out[k2] = v2
        elif type(v) in [int, float, np.ndarray]:
            print("  {}".format(v))
            gfile_out[k] = v
    pickle.dump(gfile_out, f)
print("Saved to {}".format(fname))
print(gfile_out["RMAXIS"])
"""
#gfile_path ="/home/zamp/data/200083_2400.pickle"
#gfile_path ="/home/zamp/data/201661_3000.pickle"
#gfile_path ="/home/zamp/data/203188_3000.pickle"
#gfile_path ="/home/zamp/data/200682_3500.pickle"
#gfile_path ="/home/zamp/data/193765_5500.pickle"
gfile_path ="/home/zamp/data/170361_4300.pickle"
with open(gfile_path, "rb") as f:
	gfile = pickle.load(f)

# Load wall coordinates. Run the following command beforehand:
#  sudo mount -t drvfs G: /mnt/g
xlpath = "/mnt/g/My Drive/Research/Documents/2025/excel_wall_erosion_v2.xlsx"
df = pd.read_excel(xlpath, sheet_name="wall_data_v2")

# Arrays for convienence
r1 = df["R1"].values
z1 = df["Z1"].values
r2 = df["R2"].values
z2 = df["Z2"].values
rmid = df["mid_R"].values
zmid = df["mid_Z"].values
region = df["region"].values
R, Z = np.meshgrid(gfile["R"], gfile["Z"])
Bp = gfile["Bp"]
Bz = gfile["Bz"]
Br = gfile["Br"]
Bt = gfile["Bt"]
psin = gfile["PSIRZ_NORM"]
rlim = gfile["RLIM"]
zlim = gfile["ZLIM"]
B = np.sqrt(np.square(Br) + np.square(Bz) + np.square(Bt))

# Find nearest separatrix value at the OMP
sepr = gfile["RBBBS"]
sepz = gfile["ZBBBS"]
rmaxis = gfile["RMAXIS"]
zmaxis = gfile["ZMAXIS"]

# Calculate distance to Z position of magentic axis and get sorted indices
zdist = np.abs(sepz - zmaxis)
sort_idx = np.argsort(zdist)

# Grab the nearest value where sepr > rmaxis
for idx in sort_idx:
	if sepr[idx] > rmaxis:
		rsepomp = sepr[idx]
		zsepomp = sepz[idx]
		break

print("Separatrix OMP Values")
print("  Rsep = {:.3f}".format(rsepomp))
print("  Zsep = {:.3f}".format(zsepomp))

# Now find the components of the magnetic field at that location
dist = np.sqrt(np.square(R - rsepomp) + np.square(Z - zsepomp))
min_index = np.unravel_index(np.argmin(dist), dist.shape)
#print("  Bp = {:.3f}".format(Bp[min_index]))
#print("  Br = {:.3f}".format(Br[min_index]))
#print("  Bt = {:.3f}".format(Bt[min_index]))

wall_bp = []
wall_b = []
wall_psin = []
cos_theta = []
sin_theta = []
impact_ang = []
for i in range(0, len(rmid)):
	
	# Calculate distance to each gfile data point and get index of minimum
	dist = np.sqrt(np.square(R - rmid[i]) + np.square(Z - zmid[i]))
	min_index = np.unravel_index(np.argmin(dist), dist.shape)

	# Get Bp and B
	wall_bp.append(Bp[min_index])
	wall_b.append(B[min_index])

	# Get psin
	wall_psin.append(psin[min_index])

	# Calculate components of normal vector to surface (line segment)
	Nr = z2[i] - z1[i]
	Nz = r1[i] - r2[i]  # Note sign difference
	Nt = 0.0            # Assume toroidally symmetric for now
	N = np.sqrt(np.square(Nr) + np.square(Nz) + np.square(Nt))

	# Calculate cos(theta), angle between normal vector and magnetic field
	ct = (Br[min_index] * Nr + Bz[min_index] * Nz) / (B[min_index] * N)
	cos_theta.append(ct) 

	# What we actually want though is the angle between the surface and the 
	# magnetic field, which is sin(theta)
	ang = np.arccos(ct)
	sin_theta.append(np.sin(ang))
	impact_ang.append(np.pi/2 - ang)
	

# Calculate what R-Rsep @ OMP is for each psin value at the wall. To do this,
# we first find what z index in the 2D data corresponds to the OMP (i.e., what
# (row of cells in the 2D plots).
omp_zidx = np.argmin(np.abs(gfile["Z"] - zmaxis))

# We create a mask to remove the inboard data since we only want to map to 
# the OMP.
omp_mask = gfile["R"] > rmaxis

# Then index the psin data along that row. For each wall segment, find the 
# nearest psin_omp value and corresponding Romp value.
psin_omp = psin[omp_zidx][omp_mask]
romp = gfile["R"][omp_mask]

wall_rmrsomp = []
for p in wall_psin:
	idx = np.argmin(np.abs(p - psin_omp))
	wall_romp = romp[idx]
	wall_rmrsomp.append(wall_romp - rsepomp)

# Just list what you can print out on your own
print("Data to print out")
print("  wall_bp")
print("  wall_b")
print("  wall_psin")
print("  cos_theta")
print("  sin_theta")
print("  impact_ang")
print("  wall_rmrsomp")
print("  sepr")
print("  sepz")
