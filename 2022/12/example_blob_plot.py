import pandas as pd
import matplotlib.pyplot as plt


path = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190486_1.tab"
df = pd.read_csv(path, delimiter="\t")

rmrs = df["R-Rsep(cm)"]
vr = df["Vr(m/s)"]
timestep = (df["Time(ms)"].iloc[1] - df["Time(ms)"].iloc[0]) / 1000
fblob = df["Npeaks"] / timestep
dt = df["T_blob(e-6s)"] * 1e-6
vr_avg = vr * fblob * dt
drad = df["D_rad(cm)"]

lw = 2
color = "tab:red"
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(5, 6), sharex=True)
ax1.plot(rmrs, vr, lw=lw, color=color)
ax2.plot(rmrs, dt, lw=lw, color=color)
ax3.plot(rmrs, fblob, lw=lw, color=color)
ax4.plot(rmrs, vr_avg, lw=lw, color=color)
ax1.set_ylabel(r"$\mathdefault{v_r\ (m/s)}$")
ax2.set_ylabel(r"$\mathdefault{\Delta t\ (s)}$")
ax3.set_ylabel(r"$\mathdefault{f_{blob}\ (Hz)}$")
ax4.set_ylabel(r"$\mathdefault{\overline{v_r}\ (m/s)}$")
ax4.set_xlabel(r"$\mathdefault{R-R_{sep}\ (cm)}$")
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
fig.tight_layout()
fig.show()

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(rmrs, drad, lw=lw, color=color)
ax.set_ylabel("Blob Width (cm)")
ax.grid()
ax.set_xlabel(r"$\mathdefault{R-R_{sep}\ (cm)}$")
fig.tight_layout()
fig.show()
