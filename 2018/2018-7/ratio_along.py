import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# All the H-mode, forward Bt, single null probes with enough content on them.
plist = ['A15', 'A17', 'A19', 'A21', 'A28', 'A32', 'A33', 'A34', 'A35']
ratio_df = pd.DataFrame(np.zeros(100))
for p in plist:

    filename = '/home/shawn/d3dscripts/Data/CP Excel Sheets/' + p + '.xlsx'
    df = pd.read_excel(filename)

    r_ompD = df['R-Rsep omp D (cm)'].values
    r_ompU = df['R-Rsep omp U (cm)'].values
    wD = df['W Areal Density D (1e15 W/cm2)'].values
    wU = df['W Areal Density U (1e15 W/cm2)'].values

    # Replace zeros with a really small number to prevent errors.
    #small_num = 0.00000000001
    #wD[np.where(wD==0)] = small_num
    #wU[np.where(wU==0)] = small_num

    # Create an interpolation function between the max of the min of the two,
    # and the min of the max of the two.
    fD = interp1d(r_ompD, wD)
    fU = interp1d(r_ompU, wU)

    common_x = np.linspace(np.max((r_ompD.min(), r_ompU.min())), np.min((r_ompD.max(), r_ompU.max())), 100)

    # Get the ITF/OTF ratio. For forward Bt, this is U/D. For reverse it's D/U.
    # A31 and 32 are the reverse H-Mode ones.
    ratio = fU(common_x) / fD(common_x)

    tmp_ratio = ratio
    ratio = ratio[np.logical_not(np.isnan(tmp_ratio))]
    common_x = common_x[np.logical_not(np.isnan(tmp_ratio))]
    tmp_ratio = ratio
    ratio = ratio[np.logical_not(np.isinf(tmp_ratio))]
    common_x = common_x[np.logical_not(np.isinf(tmp_ratio))]


    # Put ratios into df.
    ratio_df[p + ' R-Rsep omp (cm)'] = pd.Series(common_x)
    ratio_df[p + ' Ratio'] = pd.Series(ratio)

ratio_df.to_excel('ratios_along.xlsx')
