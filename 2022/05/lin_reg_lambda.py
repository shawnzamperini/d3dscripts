# Do a linear regression on the lamnda_ne data.
from sklearn.linear_model import LinearRegression
import pandas as pd
import matplotlib.pyplot as plt


xlpath = "/Users/zamperini/My Drive/Research/Documents/2022/05/tile_freqs.xlsx"
df = pd.read_excel(xlpath, nrows=11)

# 172413 excluded bc shoulder possibly.
shots = [167196,170843,176508,179900,164262,162766,171432,172409,174212,174305]

# Reduce DataFrame down to only shots we want.
mask = [int(s) in shots for s in df["Shot"].values]
shotdf = df[mask]
rev = [True if d=="Reverse" else False for d in shotdf["BT"]]

X = shotdf[["IP", "Peak Frequency"]].values
y = shotdf["Lambda_ne"].values
model = LinearRegression().fit(X, y)
preds = model.predict(X)


fig, ax1 = plt.subplots()

for i in range(0, len(shotdf)):
    if rev[i]:
        ax1.scatter(preds[i], y[i], color="tab:purple")
    else:
        ax1.scatter(preds[i], y[i], color="tab:red")
ax1.plot([0, 1], [0, 1], color="k", linestyle="--")
ax1.set_xlim([0.0, 0.1])
ax1.set_ylim([0.0, 0.1])

fig.tight_layout()
fig.show()
