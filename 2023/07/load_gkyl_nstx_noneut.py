import postgkyl as pg
import numpy as np
import matplotlib.pyplot as plt


z_slice = 4
def func_data(ionDensityData):
    ionDensityInterp = pg.data.GInterpModal(ionDensityData, 1, 'ms')
    interpGrid, ionDensityValues = ionDensityInterp.interpolate()

    # get cell center coordinates
    CCC = []
    for j in range(0, len(interpGrid)):
        CCC.append((interpGrid[j][1:] + interpGrid[j][:-1]) / 2)

    x_vals = CCC[0]
    y_vals = CCC[1]
    z_vals = CCC[2]
    X, Y = np.meshgrid(x_vals, y_vals)
    ionDensityGrid = np.transpose(ionDensityValues[:, :, z_slice, 0])
    return x_vals, y_vals, X, Y, ionDensityGrid


ionDensity = "/Users/zamperini/gkyldir/nstx-noNeut-izSrc/nstx-noNeut-izSrc_electron_0.bp"
ionDensityData = pg.data.GData(ionDensity)

x_vals, y_vals, X, Y, ionDensityGrid = func_data(ionDensityData)
Z = ionDensityGrid[0,0,:,:]

fig, ax = plt.subplots()

ax.pcolormesh(X, Y, Z, shading="nearest")

fig.tight_layout()
fig.show()