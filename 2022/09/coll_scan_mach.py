# A plot looking at if the lambda_ne (representative of the radial blob transport)
# trends with the Mach number but only in higher density (collisional) conditions.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


lambdas = {
"190440_1":0.01092,
"190440_2":0.05127,
"190442_1":0.08170,
"190442_2":0.10178,
"190484_1":0.00568,
"190484_2":0.00918,
"190485_1":0.17578,
"190485_2":0.03473,
"190486_1":0.02618,
"190486_2":0.09330
}

nes = {
"190440_1":0.01092,
"190440_2":0.05127,
"190442_1":0.08170,
"190442_2":0.10178,
"190484_1":0.00568,
"190484_2":999,
"190485_1":0.17578,
"190485_2":0.03473,
"190486_1":0.02618,
"190486_2":999
}

shots = [190440, 190442, 190484, 190485, 190486]
dfs = {}
for shot in shots:
    tmp = []
    for plunge in [1, 2]:
        label = "{}_{}".format(shot, plunge)
        path = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP{}_{}.tab".format(shot, plunge)
        df = pd.read_csv(path, sep="\t")
        dfs[label] = df


# For each plunge, calculate the average Mach number between 2.29 and 2.33.
# Assume data beyond that is unreliable.
window = 2.29
end_range = 2.33
avg_m = {}
for key, df in dfs.items():
    r = df["R(cm)"].values / 100
    idx = np.where(np.logical_and(r>window, r<end_range))
    mach = df["Machn"].values[idx]
    print(key)
    print(mach)
    print()
    avg_m[key] = mach.mean()


x = list(avg_m.values())
y = list(lambdas.values())
labels = list(lambdas.keys())

fig, ax1 = plt.subplots()
ax1.scatter(x, y)
for i in range(0, len(x)):
    ax1.annotate(labels[i], (x[i], y[i]))
ax1.set_xlabel("Avg. Mach")
ax1.set_ylabel("Lambda_ne")
fig.tight_layout()
fig.show()
