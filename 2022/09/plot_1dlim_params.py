import pandas as pd
import matplotlib.pyplot as plt


path = "/Users/zamperini/My Drive/Research/Documents/2022/08/tausink.xlsx"
df = pd.read_excel(path)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13, 4))

ax1.plot(df["R_in"]*100, df["Te_in"], lw=3, color="tab:red")
ax2.plot(df["R_in"]*100, df["ne_in"], lw=3, color="tab:red")
ax3.plot(df["r"]*100, df["tausink.1"], lw=3, color="tab:red")
ax1.set_title("Te (eV)", fontsize=14)
ax2.set_title("ne (m-3)", fontsize=14)
ax3.set_title(r"$\mathdefault{\tau_{sink}\ (s)}$", fontsize=14)

fig.supxlabel("r (cm)")
fig.tight_layout()
fig.show()
