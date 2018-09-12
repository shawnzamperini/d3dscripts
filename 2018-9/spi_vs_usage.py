import pandas as pd
import wikipedia as wp
import numpy as np
import matplotlib.pyplot as plt


print("Loading SPI table...")
spi_html  = wp.page('Social Progress Index').html().encode("UTF-8")
print("Loading TPES table...")
tpes_html = wp.page('List of countries by energy consumption per capita').html().encode("UTF-8")

spi_df = pd.read_html(spi_html)[1]
spi_df.set_index(0, inplace=True)
spi_df.rename(columns=spi_df.iloc[0], inplace=True)
spi_df.drop(spi_df.index[0], inplace=True)
spi_df.dropna(inplace=True)
spi_df.drop(spi_df.index[-1])

def clip_rank(clipme):
    clipme = str(clipme)[-5:]
    return clipme

spi_df = spi_df.applymap(clip_rank)
spi_df.drop(spi_df.index[-1], inplace=True)

tpes_df = pd.read_html(tpes_html)[0]
tpes_df.set_index(0, inplace=True)
tpes_df.rename(columns=tpes_df.iloc[1], inplace=True)
tpes_df.drop(tpes_df.index[[0,1]], inplace=True)
tpes_df = tpes_df.iloc[:, [3,4,5]]  # Get only the 2013 data.

full_df = tpes_df.join(spi_df)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(full_df['kgoe/a'], full_df['2017'], '.')
ax1.set_xlim([0,5000])
fig.show()
