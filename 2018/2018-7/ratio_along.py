import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# All the H-mode, forward Bt, single null probes with enough content on them.
plist = ['A2', 'A3', 'A4', 'A7', 'A11', 'A12', 'A15', 'A17', 'A18', 'A19',
         'A20', 'A21', 'A28', 'A32', 'A33', 'A34', 'A35']
ratio_df = pd.DataFrame(np.zeros(100))
for p in plist:
#for p in ['A2']:
    try:
        print('Probe: ' + p)
        filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/' + p + '.xlsx'
        #filename = '/home/shawn/d3dscripts/Data/CP Excel Sheets/' + p + '.xlsx'
        df = pd.read_excel(filename)
        df = df[~np.isnan(df.index)]

        r_ompD = df['R-Rsep omp D (cm)'].values
        r_ompU = df['R-Rsep omp U (cm)'].values
        wD = df['W Areal Density D (1e15 W/cm2)'].values
        wU = df['W Areal Density U (1e15 W/cm2)'].values

        if p == 'A2':
            r_ompU = r_ompU[:-1]
            wU = wU[:-1]

        # Filter out low signal points.
        thresh = 0.002
        tmp_wD = np.copy(wD)
        tmp_wU = np.copy(wU)
        wD = tmp_wD[tmp_wD > thresh]
        wU = tmp_wU[tmp_wU > thresh]
        r_ompD = r_ompD[tmp_wD > thresh]
        r_ompU = r_ompU[tmp_wU > thresh]

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
        if p in ['A2', 'A3', 'A4', 'A7', 'A8', 'A31', 'A32']:
            print("Ratio is D/U")
            ratio = fD(common_x) / fU(common_x)
        else:
            print("Ratio is U/D")
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

    except:
        print("Error on probe " + p)

ratio_df.to_excel('ratios_along.xlsx')
