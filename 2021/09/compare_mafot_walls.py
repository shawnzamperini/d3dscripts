import matplotlib.pyplot as plt
import pandas as pd

col_names = ["R[m]", "Z[m]", "N_toroidal", "connection length [km]",
  "psimin (penetration depth)", "psimax", "psiav", "FL pitch angle",
  "FL yaw angle", "theta", "psi"]
root = "/Users/zamperini/Documents/d3d_work/"

def read_df(fname):
    df = pd.read_csv(root+fname, delimiter="\t", skiprows=52, names=col_names)
    return df

wall_2d = read_df("lam_184527_both_2d_wall.dat")
wall_3d = read_df("lam_184527_both_3d_wall.dat")
wall_2d_p1 = read_df("lam_184527_+1_2d_wall.dat")
wall_2d_m1 = read_df("lam_184527_-1_2d_wall.dat")
wall_3d_p1 = read_df("lam_184527_+1_3d_wall.dat")
wall_3d_m1 = read_df("lam_184527_-1_3d_wall.dat")

df = read_df("lam_187122_+1_3d_wall.dat")

fig, ax = plt.subplots()
ax.plot(wall_2d["R[m]"], wall_2d["connection length [km]"]*1000, "tab:red", label="2D")
ax.plot(wall_3d["R[m]"], wall_3d["connection length [km]"]*1000, "tab:purple", label="3D")
#ax.plot(wall_2d_p1["R[m]"], wall_2d_p1["connection length [km]"]*1000, "tab:red", label="2D OTF", linestyle="-", lw=2)
#ax.plot(wall_2d_m1["R[m]"], wall_2d_m1["connection length [km]"]*1000, "tab:red", label="2D ITF", linestyle="-.", lw=2)
#ax.plot(wall_3d_p1["R[m]"], wall_3d_p1["connection length [km]"]*1000, "tab:purple", label="3D OTF", linestyle="-", lw=2)
#ax.plot(wall_3d_m1["R[m]"], wall_3d_m1["connection length [km]"]*1000, "tab:purple", label="3D ITF", linestyle="-.", lw=2)
#ax.plot(df["R[m]"], df["connection length [km]"]*1000, "k")
ax.set_xlabel("R (m)", fontsize=14)
ax.set_ylabel("Connection Length (m)", fontsize=14)
ax.grid()
ax.set_yscale("log")
ax.legend(fontsize=14)
fig.tight_layout()
fig.show()
