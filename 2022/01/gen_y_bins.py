import numpy as np


ystart = 0.0025
yend = 44.00
ny = 150
first_step = 0.0025
last_step = 1.00
num_ys = 150

mults = np.linspace(1e-10, 1.0, 10000)

# generate bins with the starting multiplier for the step increase.
prev_vals = []
vals = []
for i in range(0, len(mults)):
    prev_vals = np.copy(vals)
    vals = [ystart]
    y = vals[0]
    j = 1
    step = first_step
    while y < yend:
        vals.append(vals[j-1] + step)
        y = vals[-1]
        step += mults[i] * step
        j += 1
    if len(vals) > num_ys:
        continue
    else:
        break

for k in prev_vals:
    print("{:.4}".format(k))
