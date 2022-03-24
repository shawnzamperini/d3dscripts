import numpy as np
import matplotlib.pyplot as plt
import random
from tqdm import tqdm


nsims = 1000
batch1_fracs = np.linspace(0.01, 0.99, 100)
sims_needed = np.zeros(len(batch1_fracs))

for i in tqdm(range(0, len(batch1_fracs))):

    for j in range(0, nsims):
        count = 0
        napples = 100
        while napples > 1:

            # Probability batch 1 has the bad apple is just the fraction.
            bad_apple = random.uniform(0, 1)
            if bad_apple < batch1_fracs[i]:

                # If batch1 had the bad apple, then remove then good apples from the
                # pool and restart.
                napples = napples - napples * (1 - batch1_fracs[i])

            # If batch 2 had the bad apple, remove all of batch 1 from the pool.
            else:
                napples = napples - napples * batch1_fracs[i]
            count += 1

        sims_needed[i] += count / nsims

best_x = batch1_fracs[np.argmin(sims_needed)]
best_y = np.min(sims_needed)

plt.xkcd()
fig, ax = plt.subplots()
ax.plot(batch1_fracs, sims_needed, lw=3, color="k")
ax.annotate("Best fraction = {:.2f}, {:d} splits".format(best_x, int(best_y)),
  (best_x, best_y), (0.3, 0.5), textcoords="figure fraction",
  arrowprops=dict(facecolor='black', shrink=0.05, width=2))
ax.set_xlabel("Fraction of apples in batch 1")
ax.set_title("Average number of splits until\nbad apple identified")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
fig.tight_layout()
fig.show()
