# Some work to looks for trends in some of Jessica's datasets.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split


path1 = "/Users/zamperini/Downloads/mortality over time with average island rainfall.csv"
df1 = pd.read_csv(path1)


# One of the big dependent variables is survivorship.
Xcols = ["rainfall", "Elevation_m", "Slope_degrees", "days.since.sow"]
X = df1.dropna()[Xcols]
y = df1.dropna()["survivorship"]
X_train, X_test, y_train, y_test = train_test_split(X, y)

# Standardize the data...

# Train and score...
model = LinearRegression().fit(X_train, y_train)


y_pred = model.predict(X_test)

fig, ax1 = plt.subplots()
ax1.scatter(y_pred, y_test)
ax1.set_xlabel("predicted")
ax1.set_ylabel("actual")
fig.tight_layout()
fig.show()
