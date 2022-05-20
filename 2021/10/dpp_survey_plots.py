import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams['font.family'] = 'Century Gothic'

# Load, select just columns with numbers, set question number as column.
#path = "/Users/zamperini/Documents/d3d_work/211025 DEI Survey (Responses).xlsx"
path = "/Users/zamperini/Documents/d3d_work/dei/211025 DEI Survey (Responses).xlsx"
df = pd.read_excel(path)

answer_cols = df.columns[[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 24]]
answer_nums = np.arange(1, 14, dtype=int)

# Instead choose a representative subset.
#answer_cols = df.columns[[1, 7, 13, 17]]
#answer_nums = np.arange(1, 5, dtype=int)

answers = df[answer_cols]
answers.columns = answer_nums

# Instead choose a representative subset.
numeric_answers = np.arange(1, 11, dtype=int)
#numeric_answers = np.array([1, 4, 7, 8])

min_answers = answers[answers[13] == "Yes"][numeric_answers]
maj_answers = answers[answers[13] == "No"][numeric_answers]


# Average repsonses difference.
avg_diff = ((maj_answers.mean() - min_answers.mean()) / maj_answers.mean()).mean()
print("Minorities responded on average {:.1f}% lower".format(avg_diff * 100))

xlabels = ["Diversity is valued", "Diversity is valued\nby leadership",
  "Diversity is valued\nby supervisor", "Respected by\ncolleagues",
  "Identity is\nvalued", "Equal\nopportunity", "Supported in\ncareer",
  "Discussing DEI", "Good place\nto work", "Trust in\nhandling issues"]
#xlabels = ["Diversity is valued\nat DIII-D", "I am respected by\nmy colleagues",
#  "Everyone is\nsupported in their\ncareers at DIII-D",
#  "I feel comfortable\ndiscussing DEI issues\nwith my colleagues"]
xlabels = np.array(xlabels)


# Bar graph of the answers separated by minority or not.
min_avg = min_answers.mean().values
maj_avg = maj_answers.mean().values
sort_idx = np.argsort(min_avg)[::-1]
x = np.arange(len(min_avg))
width = 0.35
fig, ax = plt.subplots(figsize=(9, 5))
ax.grid(zorder=0)
ax.bar(x - width/2, maj_avg[sort_idx], width, label="Well-represented", color="tab:red", zorder=50)
ax.bar(x + width/2, min_avg[sort_idx], width, label="Under-represented", color="tab:purple", zorder=51)
ax.set_xticks(x)
#ax.set_xticklabels(xlabels, rotation=0, ma="left", fontsize=14)
ax.set_xticklabels(xlabels[sort_idx], rotation=270, ha="center", fontsize=14)
ax.legend(fontsize=14, framealpha=1.0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim([1, 5])
ax.set_yticks(np.arange(1, 6, dtype=int))
ylabels = ["Strongly\nDisagree", "Disagree", "Neutral", "Agree", "Strongly\nAgree"]
ax.set_yticklabels(ylabels, fontsize=14)
ax.set_ylabel("Mean Response", fontsize=14)
#ax.axhline(1.0, linestyle="--", color="k")
#ax.axhline(3.0, linestyle="--", color="k")
#ax.axhline(5.0, linestyle="--", color="k")

fig.tight_layout()
fig.show()

# Plot distribution of a given question.
min_dist = min_answers[1]
maj_dist = maj_answers[1]
fig, ax = plt.subplots()
maj_bins = ax.hist(maj_dist, [1, 2, 3, 4, 5, 6], edgecolor="tab:red", label="Well-represented", density=True, fc=(1, 0, 0, 0), rwidth=0.75)
min_bins = ax.hist(min_dist, [1, 2, 3, 4, 5, 6], edgecolor="tab:purple", label="Under-represented", density=True, fc=(1, 0, 0, 0), rwidth=0.75)
ax.legend(fontsize=14)
ax.set_xticks(np.arange(1.5, 6.5))
ax.set_xticklabels(np.arange(1, 6))
ax.set_xlabel("Response (1-5)", fontsize=14)
ax.set_ylabel("# of Responses")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()

# Pie graph showing the percentage of responders that had at least X amount
# of 1 responses.
ans = np.arange(1, 6, dtype=int)
labels = ["Strongly\nDisagree", "Disagree", "Neutral", "Agree", "Strongly\nAgree"]
min_counts = {1:0, 2:0, 3:0, 4:0, 5:0}
maj_counts = {1:0, 2:0, 3:0, 4:0, 5:0}
min_total = 0
maj_total = 0
for a in ans:
    for n in numeric_answers:
        try:
            min_counts[a] += min_answers[n].value_counts()[a]
            min_total += min_answers[n].value_counts()[a]
        except KeyError:
            pass
        try:
            maj_counts[a] += maj_answers[n].value_counts()[a]
            maj_total = maj_answers[n].value_counts()[a]
        except KeyError:
            pass
cmap = plt.get_cmap('RdYlGn')
colors = cmap(np.linspace(0, 1.0, 5))
min_x = np.array(list(min_counts.values())) / min_total
maj_x = np.array(list(maj_counts.values())) / maj_total
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
p2 = ax2.pie(min_x, labels=labels, normalize=True, textprops={"fontsize":14}, colors=colors, labeldistance=1.15)
p1 = ax1.pie(maj_x, labels=labels, normalize=True, textprops={"fontsize":14}, colors=colors, labeldistance=1.15)

# Emphasize the good and bad wedges.
#for w in range(1, 4):
#    p1[0][w].set_alpha(0.4)
#    p2[0][w].set_alpha(0.4)

ax2.set_title("Under-represented", fontsize=16, color="tab:red")
ax1.set_title("Well-represented", fontsize=16, color="tab:purple")
fig.tight_layout()
fig.show()


# Copied from https://stackoverflow.com/questions/20549016/explode-multiple-slices-of-pie-together-in-matplotlib
# redraw plot using patches (right)
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
fig, ax = plt.subplots(1, 2)
ax[0].set_aspect('equal')
ax[0].axis('off')
ax[1].set_aspect('equal')
ax[1].axis('off')

tmp_p2 = ax[1].pie(min_x, labels=labels, normalize=True, textprops={"fontsize":14}, colors=colors, labeldistance=1.25)
tmp_p1 = ax[0].pie(maj_x, labels=labels, normalize=True, textprops={"fontsize":14}, colors=colors, labeldistance=1.25)

for w in range(0, 5):
    tmp_p1[0][w].set_alpha(0.0)
    tmp_p2[0][w].set_alpha(0.0)

wedges1 = p1[0]
wedges2 = p2[0]

groups=[[0,1,2],[3,4]]
radfraction = 0.1
patches1 = []
patches2 = []
for i in groups:
  ang1 = np.deg2rad((wedges1[i[-1]].theta2 + wedges1[i[0]].theta1)/2,)
  ang2 = np.deg2rad((wedges2[i[-1]].theta2 + wedges2[i[0]].theta1)/2,)
  for j in i:
    we1 = wedges1[j]
    we2 = wedges2[j]
    center1 = (radfraction*we1.r*np.cos(ang1), radfraction*we1.r*np.sin(ang1))
    center2 = (radfraction*we2.r*np.cos(ang2), radfraction*we2.r*np.sin(ang2))
    patches1.append(mpatches.Wedge(center1, we1.r, we1.theta1, we1.theta2))
    patches2.append(mpatches.Wedge(center2, we2.r, we2.theta1, we2.theta2))

colors = np.linspace(0, 1, len(patches1))
collection1 = PatchCollection(patches1, cmap=plt.cm.RdYlGn)
collection2 = PatchCollection(patches2, cmap=plt.cm.RdYlGn)
collection1.set_array(np.array(colors))
collection2.set_array(np.array(colors))
ax[0].add_collection(collection1)
ax[0].autoscale(True)
ax[1].add_collection(collection2)
ax[1].autoscale(True)
ax[1].set_title("Under-represented", fontsize=16, color="tab:purple", weight="bold")
ax[0].set_title("Well-represented", fontsize=16, color="tab:red", weight="bold")
fig.tight_layout()
fig.show()
