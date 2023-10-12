import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_excel("/Users/zamperini/My Drive/Research/Documents/2022/02/mixed_material.xlsx",
                   sheet_name="TiB2 ZrB2 Yields")

Eimp = df["D Eimp (eV)"]
YB = df["B Y (atoms/ion)"]
YTi = df["Ti Y (atoms/ion)"]
YZr = df["Zr Y (atoms/ion)"]

f_YB = interp1d(Eimp, YB)
f_YTi = interp1d(Eimp, YTi)
f_YZr = interp1d(Eimp, YZr)

x = np.geomspace(1, 5000, 100)
y_B = f_YB(x)
y_Ti = f_YTi(x)
y_Zr = f_YZr(x)

fig, ax = plt.subplots(figsize=(5, 4))

ax.plot(x, y_B, label="B")
ax.plot(x, y_Ti, label="Ti")
ax.plot(x, y_Zr, label="Zr")

ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()
ax.set_xlabel("E Impact (eV)")
ax.set_ylabel("Yield")
fig.tight_layout()
fig.show()