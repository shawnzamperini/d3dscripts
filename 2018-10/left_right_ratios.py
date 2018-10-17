import pandas as pd
import pretty_plots as pp
import numpy as np


#probes = ['BD05', 'BU05', 'BD07', 'BU07', 'BD09', 'BU09']
#probes = ['BD05', 'BD07', 'BD09']
#probes = ['BU05', 'BU07', 'BU09']
#probes = ['CU05', 'CU07', 'CU09']
#probes = ['AU35', 'AU34']
probes = ['AD35', 'AD34']
#probes = ['AD05', 'AD06']
ratio_df = pd.DataFrame()
option = 1
if option == 0:
    for probe in probes:
        print('Now on probe ' + probe)
        filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
                   'Polodial_Scans/'+probe+'_Map.xlsx'
        sheet_name = 'Sheet2'

        df = pd.read_excel(filename, sheet_name=sheet_name).set_index('Radial [mm]')[:-1]

        # Let's look at the ratio of W between pol = 2 / pol = 0.5.
        ratio = df[np.arange(0,2.5,0.25)].sum(axis=1) / df[np.arange(2.5,5,0.25)].sum(axis=1)
        #ratio = df[2.0] / df[0.5]
        ratio_df[probe] = ratio
elif option == 1:
    for probe in probes:
        print('Now on probe ' + probe)
        filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
                   'Polodial_Scans/'+probe+'_Map.xlsx'
        sheet_name = 'Sheet1'
        df = pd.read_excel(filename, sheet_name=sheet_name).set_index('Radial [mm]')[:-1]
        df = df.drop('Poloidal [mm]')[:202]

        ratio = df[np.arange(0,1.75,0.25)].sum(axis=1) / df[np.arange(1.75,4.5,0.25)].sum(axis=1)
        #ratio = df[1.0] / df[3.75] # For A probe at least, forward
        ratio_df[probe] = ratio


all_ratios = 1.0 / ratio_df.mean(axis=1)
fig = pp.pplot(ratio_df.index.values/10, all_ratios.values, fmt='-', yrange=[0,2.5],
               xlabel='Radial (cm)', ylabel='Top Half / Bottom Half')
