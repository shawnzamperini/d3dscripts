# See if any obvious scalings from combining variables comes out for the LP peaks.
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt


# Load in, restrict to just L-mode.
path = "/Users/zamperini/My Drive/Research/Documents/2022/09/lp_peak_scan.xlsx"
df = pd.read_excel(path, sheet_name="grouping")

ldf = df[df["mode"]=="l"]
ldf = ldf[ldf["include_python"]=="y"]
ldf = ldf[["peak location (cm)", "ip", "fg", "lambda_ne", "flux_expansion"]].dropna()
#fitvars = ["lambda_ne"]
fitvars = ["lambda_ne", "flux_expansion"]
X = ldf[fitvars]
y = ldf["peak location (cm)"]

# Standardize the data.
scaler = StandardScaler().fit(X)
X_scaled = scaler.transform(X)

if len(fitvars) == 2:
    # Try a power law.
    def power_law(xx, a, b, c):
        return a * xx[0]**b * xx[1]**c
    popt, pcov = curve_fit(power_law, X.values.T, y)
    fity1 = power_law(X.values.T, *popt)
    eq1_str = "{:.2f}fg^{:.2f} lambda_ne^{:.2f}".format(popt[0], popt[1], popt[2])

# Try a linear regression.
model = LinearRegression().fit(X_scaled, y)
fity2 = model.predict(X_scaled)
#eq2_str = "{:.2f}fg + {:.2f}lambda_ne".format(model.coef_[0], model.coef_[1])


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

confx = np.linspace(-10, 10, 3)

if len(fitvars) == 2:
    ax1.set_facecolor("tab:red")
    ax1.fill_between(confx, confx*2, confx/2, zorder=1, alpha=0.75, color="tab:green")
    ax1.plot([-10, 10], [-10, 10], color="k", linestyle="--", lw=2, zorder=5)
    ax1.scatter(fity1, y, color="k", zorder=15)
    ax1.grid(zorder=3)
    ax1.set_xlabel(eq1_str, fontsize=14)
    ax1.set_ylabel("Peak Location (cm)", fontsize=14)
    ax1.set_xlim([0, 10])
    ax1.set_ylim([0, 10])

ax2.set_facecolor("tab:red")
ax2.fill_between(confx, confx*2, confx/2, zorder=1, alpha=0.75, color="tab:green")
ax2.plot([-10, 10], [-10, 10], color="k", linestyle="--", lw=2, zorder=5)
ax2.scatter(fity2, y, color="k", zorder=15)
ax2.grid(zorder=3)
ax2.set_xlabel("Predicted Location (cm)", fontsize=14)
ax2.set_ylabel("Peak Location (cm)", fontsize=14)
ax2.set_xlim([0, 10])
ax2.set_ylim([0, 10])

fig.tight_layout()
fig.show()
