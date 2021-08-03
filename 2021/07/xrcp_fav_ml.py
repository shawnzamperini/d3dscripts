import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split


all_df = pd.read_excel("xrcp_data.xlsx")
fav_df = all_df[all_df["direction"] == "unfavorable"]


fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

color_boss = "pinj"
fav_df["color_scale"] = (fav_df[color_boss] - fav_df[color_boss].min()) / + \
  (fav_df[color_boss].max() - fav_df[color_boss].min())
cmap = plt.get_cmap('Reds')

avg_machs = []; avg_boss = []
for id in fav_df["id"].unique():
    df = fav_df[fav_df["id"] == id]
    color = cmap(df["color_scale"].mean())
    ax1.plot(df["Rho"], df["Machn"], color=color)

    mask = np.logical_and(df["Rho"] >= 1.02, df["Rho"] <= 1.05)
    avg_machs.append(df["Machn"][mask].mean())
    avg_boss.append(df[color_boss].mean())

ax1.set_ylim([-1, 1])
ax1.set_xlabel("Rho")
ax1.set_ylabel("Mach")

fig.tight_layout()
fig.show()

clean = np.logical_and(fav_df["Rho"] >= 1.02, fav_df["Rho"] <= 1.05)
clean_df = fav_df[clean]
avg_df = pd.DataFrame()
for id in clean_df["id"].unique():
    avg_df = avg_df.append(clean_df[clean_df["id"] == id].mean(), ignore_index=True)

X = avg_df[['tinj', "densv2", 'Ne (E18 m-3)', 'Te(eV)', 'ne_ped', 'pinj', 'psol', 'te_ped']]
y = avg_df["Machn"]
X_train, X_test, y_train, y_test = train_test_split(X, y)

scaler = StandardScaler().fit(X_train)
X_train_scaled = scaler.transform(X_train)
X_test_scaled = scaler.transform(X_test)

model = LinearRegression().fit(X_train_scaled, y_train)
print("Train score: {:.3f}".format(model.score(X_train_scaled, y_train)))
print("Test score:  {:.3f}".format(model.score(X_test_scaled, y_test)))
