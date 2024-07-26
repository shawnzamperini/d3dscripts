import pandas as pd
import matplotlib.pyplot as plt


path = "/Users/zamperini/My Drive/Research/Documents/2022/10/sol29_vrs.xlsx"
df = pd.read_excel(path)

fig, ax = plt.subplots(figsize=(5,4))

ax.hist(df["vr"], bins=50, density=True, label="SOL 29")
ax.plot(df["vr-scipy"], df["prob"], label="scipy")
ax.set_xlabel("vr (m/s)")
ax.set_ylabel("probability")
ax.legend()

fig.tight_layout()
fig.show()
