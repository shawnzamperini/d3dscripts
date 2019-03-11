import pretty_plots as pp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


df = pd.read_csv('grid.txt', sep=' ', index_col=0)
df.drop(df.columns[-1], axis=1, inplace=True)

x = np.array(df.index.values,   dtype=np.float64)
y = np.array(df.columns.values, dtype=np.float64)

X, Y = np.meshgrid(x, y)
Z = df.values.T

pp.ppcontourf(Y, X, Z, xlabel='Parallel (m)', ylabel='Radial (m)', cbarlabel='ne (m-3)', yrange=[x.max(), x.min()], cmap='viridis')

#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.contourf(X, Y, Z)
#fig.show()
