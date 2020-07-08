# Needed libraries.
import pandas as pd
from scipy.signal import find_peaks


# Constants telling us the threshold values for peak detection.
site_thresh = {'GC':0.4, 'CLM':0.4, 'SLM':0.4, 'GATR':0.6, 'BOXTREE':0.22, 'WF':0.18, 'WF0':0.18, 'INDIANTOWN':0.45, 'UPC':0.3, 'OYSTER':0.23}

# Load in DataFrame.
water_path = '/mnt/c/Users/Shawn/Documents/water_level_data.xlsx'
print("Loading {}...".format(water_path))
water_df = pd.read_excel(water_path)

# A DataFrame that will hold all of our average values.
avg_df = pd.DataFrame(columns=['date', 'site', 'depth', 'changeavg_df'])

# Go through one site at a time.
for site in water_df['site'].unique():

    # Printout.
    print("Site: {}".format(site))

    # Use subset of the DataFrame for just this site.
    site_df = water_df[water_df['site'] == site]

    # Find all the peaks. Require that they be greater than a site-specific
    # threshold and be at least 10 data points apart.
    idx = find_peaks(site_df['depth'], height=site_thresh[site], distance=10)

    # Go through one date at a time and find the average peak heights.
    for date in site_df['time'].dt.date.unique():
        avg_depth = site_df[site_df['time'].dt.date == date]['depth'].mean()
        max_depth = site_df[site_df['time'].dt.date == date]['depth'].max()
        min_depth = site_df[site_df['time'].dt.date == date]['depth'].min()
        depth_change = max_depth - min_depth

        # Store the results in our avg_df.
        avg_df = avg_df.append({'date':date, 'site':site, 'depth':avg_depth, 'change':depth_change}, ignore_index=True)

# Save to an Excel file.
avg_df.to_excel('avg_water.xlsx')
