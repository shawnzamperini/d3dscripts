import numpy as np
import matplotlib.pyplot as plt


near_x  = np.linspace(4, 11.5, 100)
near_te = 25.185 * np.exp(-0.252 * near_x)
near_ne = 22.602 * np.exp(-0.239 * near_x)
far_x   = np.linspace(11.5, 14, 50)
far_te  = 23880  * np.exp(-0.845 * far_x)
far_ne  = 136568 * np.exp(-0.979 * far_x)

# Force them to match.
far_te = far_te * near_te[-1] / far_te[0]
far_ne = far_ne * near_ne[-1] / far_ne[0]

x = np.append(near_x, far_x)
te = np.append(near_te, far_te)
ne = np.append(near_ne, far_ne) * 10**(18)

te[te < 5.0] = 5.0

# 2L ~ 10 from 5- 8.5 cm, then goes down to ~3 from 8.5 to 11, then ~0.5
# outwards. These are pretty rough approximations.
l = np.zeros(len(x))
for i in range(0, len(l)):
    if x[i] < 11.5:
        l[i] = 10.0 / 2.0
    else:
        l[i] = 3.0 / 2

nu = 10**(-16) * 5.0 * ne / (te**2)
#nu = 10**(-16) * l * ne / (te**2)

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

fig, host = plt.subplots(figsize=(8, 5))
fig.subplots_adjust(right=0.75)

par1 = host.twinx()
par2 = host.twinx()

# Offset the right spine of par2.  The ticks and label have already been
# placed on the right by twinx above.
par2.spines["right"].set_position(("axes", 1.2))
# Having been created by twinx, par2 has its frame off, so the line of its
# detached spine is invisible.  First, activate the frame but make the patch
# and spines invisible.
make_patch_spines_invisible(par2)
# Second, show the right spine.
par2.spines["right"].set_visible(True)

p1, = host.plot(x, te, label="Te", color="r", lw=4)
p2, = par1.plot(x, ne, label="ne", color="b", lw=4)
p3, = par2.plot(x, nu, label=r"$\mathdefault{\nu_{SOL}^*}$", color="g", lw=4)

host.set_xlabel("R-Rsep OMP (cm)", fontsize=16)
host.set_ylabel("Te (eV)", fontsize=16)
par1.set_ylabel("ne (m-3)", fontsize=16)
par2.set_ylabel(r"$\mathdefault{\nu_{SOL}^*}$", fontsize=16)

#host.set_yscale("log")
#par1.set_yscale("log")
#par2.set_yscale("log")

host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())
par2.yaxis.label.set_color(p3.get_color())

tkw = dict(size=4, width=1.5)
host.tick_params(axis='y', which="both", colors=p1.get_color(), **tkw, labelsize=12)
par1.tick_params(axis='y', which="both", colors=p2.get_color(), **tkw, labelsize=12)
par2.tick_params(axis='y', which="both", colors=p3.get_color(), **tkw, labelsize=12)
host.tick_params(axis='x', which="both", **tkw, labelsize=10)

host.grid()

lines = [p1, p2, p3]

host.legend(lines, [l.get_label() for l in lines], fontsize=12)
fig.tight_layout()
fig.show()

#fig, (ax1, ax2) = plt.subplots(1, 2)
#ax1.plot(x, te, label="Te")
#ax1b = ax1.twinx()
#ax1b.plot(x, ne, label="ne")
#ax1.legend()
#x2.plot(x, nu)
#ax2.set_xlabel("R-Rsep OMP (cm)", fontsize=16)
#x2.set_ylabel("Collisionality", fontsize=16)
#fig.tight_layout()
#fig.show()
