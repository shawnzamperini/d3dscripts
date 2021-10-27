# This script is part of the DIII-D DEI efforts where the DEI panel filled out
# the NERCHE rubric adapted for DIII-D. Some plots of the data are made here.
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


results_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/NERCHE_responses.xlsx"
results = pd.read_excel(results_path, sheet_name="Data", usecols="B:AC", nrows=15)
nums = np.arange(1, 29, dtype=int)

# Bar graph of the mean results sorted.
mean_vals = results.mean().values
sort_idx = np.argsort(mean_vals)
fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(nums, mean_vals[sort_idx], zorder=3, edgecolor="k", color="tab:red")
ax.set_xticks(nums)
ax.set_xticklabels(nums[sort_idx], rotation=-45)
ax.set_xlabel("Question Number")
ax.set_ylabel("Average Stage")
ax.set_ylim(0, 3)
ax.grid(axis="y", zorder=2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()

# Print out the average results of each cetegory.
phil = results.columns[1:3]
fclty = results.columns[3:8]
teach = results.columns[8:11]
staff = results.columns[11:14]
grad = results.columns[14:17]
admin = results.columns[17:]
overall = results.columns

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
