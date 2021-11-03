# This script is part of the DIII-D DEI efforts where the DEI panel filled out
# the NERCHE rubric adapted for DIII-D. Some plots of the data are made here.
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from cycler import cycler


results_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/NERCHE_responses.xlsx"
results = pd.read_excel(results_path, sheet_name="Data", usecols="B:AC", nrows=15)
nums = np.arange(1, 29, dtype=int)

# Groups that each question is a part of.
phil = results.columns[1:3]
fclty = results.columns[3:8]
teach = results.columns[8:11]
staff = results.columns[11:14]
grad = results.columns[14:17]
admin = results.columns[17:]
overall = results.columns

cmap = plt.get_cmap('Set1')
colors = cmap(np.linspace(0, 1, 9))
bar_colors = [0, 1, 2, 3, 4, 7]

fig, ax = plt.subplots(figsize=(10, 4))
prev_x = 0
count = 0
for cols in [phil, fclty, teach, staff, grad, admin]:
    x = np.arange(0, len(cols)) + prev_x + 2
    mean_vals = results[cols].mean().values
    sort_idx = np.argsort(mean_vals)
    prev_x = x[-1]
    ax.bar(x, mean_vals[sort_idx], zorder=3, edgecolor="k", color=colors[bar_colors[count]])
    #ax.set_xticks(x)
    count += 1
ax.set_ylabel("Average Stage", fontsize=14)
ax.set_ylim(0, 3)
ax.grid(axis="y", zorder=2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xticks([])
fig.tight_layout()
fig.show()

def get_mean(cols):
    sum = results[cols].sum().sum()
    count = results[cols].count().sum()
    return sum / count

print("Philosophy: {:.2f}".format(get_mean(phil)))
print("Faculty suport for DEI: {:.2f}".format(get_mean(fclty)))
print("Teaching: {:.2f}".format(get_mean(teach)))
print("Staff engagement: {:.2f}".format(get_mean(staff)))
print("Graduate/postdoc: {:.2f}".format(get_mean(grad)))
print("Administrative: {:.2f}".format(get_mean(admin)))
print("Overall: {:.2f}".format(get_mean(overall)))
