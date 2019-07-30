import numpy  as np
import pandas as pd
import pretty_plots as pp


# Paths to the Excel files for each side.
a08itf_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/AD08_Map_Analysis.xlsx"
a08otf_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/AU08_Map_Analysis.xlsx"
a15itf_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/AU15_Map_Analysis.xlsx"
a15otf_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/AD15_Map_Analysis.xlsx"

#b05itf_path  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/BD05_Map_Analysis.xlsx"
b05itf_path  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/BD05_Map_Incomplete_Analysis.xlsx"
b05otf_path  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/BU05_Map_Analysis2.xlsx"
b07itf_path  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/BU07_Map_Analysis.xlsx"
b07otf_path  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/BD07_Map_Analysis.xlsx"

c05itf_path  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/CD05_Map_Analysis2.xlsx"
c05otf_path  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/CU05_Map_Analysis2.xlsx"
c07itf_path  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/CU07_Map_Analysis.xlsx"
c07otf_path  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/CD07_Map_Analysis.xlsx"

# Load the RBS files.
a08rbs_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A8.xlsx"
a15rbs_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A15.xlsx"
b05rbs_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/B5.xlsx"
c05rbs_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/C5.xlsx"
b07rbs_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/B7.xlsx"
c07rbs_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/C7.xlsx"


def plot_with_rbs(lams_path, rbs_path, u_or_d, cal_slope, cal_intercept, middle=2.0, r_shift=0, color=6, avg=False):
    """
    u_or_d: Either 'U' or 'D'.
    """

    # Load DataFrames.
    lams_df = pd.read_excel(lams_path, sheet_name='MapData')
    rbs_df  = pd.read_excel(rbs_path)

    # Get relevant RBS data.
    dist = rbs_df['Distance from Tip ' + u_or_d + ' (cm)'].values
    omp  = rbs_df['R-Rsep omp ' + u_or_d + ' (cm)'].values
    w_rbs    = rbs_df['W Areal Density ' + u_or_d + ' (1e15 W/cm2)'].values

    # Do a fit to get R-Rsep OMP for the LAMS data.
    p = np.polyfit(dist, omp, 1)
    lams_df['R-Rsep omp (cm)'] = p[1] + p[0] * lams_df['Axial Location [mm]'] / 10.0

    # Create a 2D DataFrame of just the Total W measurements (for RBS comparison).
    omp_locs = lams_df['R-Rsep omp (cm)'].unique()
    z_locs   = lams_df['z Location [mm]'].unique()
    total_df = pd.DataFrame()
    for z in z_locs:
        total_df[z] = lams_df[lams_df['z Location [mm]'] == z]['Total W'].values

    # Set index as the omp locations.
    total_df.set_index(omp_locs, inplace=True)

    if avg:
        r = total_df[middle].index.values
        w = total_df.median(axis=1).values
    else:
        # Get the centerline data near where RBS was taken.
        r = total_df[middle].index.values
        w = total_df[middle].values

    # Use calibration to convert from counts to areal density.
    w = w * cal_slope + cal_intercept

    # Shift LAMS R data if asked.
    r = r + r_shift

    # Plot it.
    fig = pp.pplot(r, w, fmt='-', label='LAMS', color=color)
    fig = pp.pplot(omp, w_rbs, fmt='.', ms=15, fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)', label='RBS', color=color)

    return {'LAMS Romp':r, 'LAMS W':w}


def run(probe):
    """
    probe: Either 'a', 'b', or 'c'.
    """

    if probe == 'a':

        middle = 2.0

        a08itf = plot_with_rbs(a08itf_path, a08rbs_path, 'D', 0.5E-06, 0, color=8, middle=middle)
        a08otf = plot_with_rbs(a08otf_path, a08rbs_path, 'U', 0.5E-06, 0, color=8, middle=middle)
        a15itf = plot_with_rbs(a15itf_path, a15rbs_path, 'U', 5.015E-07, 0, r_shift=0, middle=middle)
        a15otf = plot_with_rbs(a15otf_path, a15rbs_path, 'D', 5.015E-07, 0, r_shift=0, middle=middle)

        ignore = 15

        fig = pp.pplot(a15itf['LAMS Romp'][ignore+55:], a15itf['LAMS W'][ignore+55:], fmt='-',  color=6, lw=7)
        fig = pp.pplot(a15otf['LAMS Romp'][ignore+90:], a15otf['LAMS W'][ignore+90:], fmt='-', color=6, lw=3, fig=fig)
        fig = pp.pplot(a08itf['LAMS Romp'][ignore:], a08itf['LAMS W'][ignore:], fmt='-',  color=8, lw=7, fig=fig)
        fig = pp.pplot(a08otf['LAMS Romp'][ignore:], a08otf['LAMS W'][ignore:], fmt='-', color=8, lw=3, fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)')

    elif probe == 'b':

        middle = 2.0

        # This is the actual calibration for B5. Very nice. Currently using same
        # value for B7 until we actually get to LAMS it.
        b05itf = plot_with_rbs(b05itf_path, b05rbs_path, 'D', 5.0E-07, 0, color=8, middle=middle)
        b05otf = plot_with_rbs(b05otf_path, b05rbs_path, 'U', 5.0E-07, 0, color=8, middle=middle)
        b07itf = plot_with_rbs(b07itf_path, b07rbs_path, 'U', 5.0E-07, 0, r_shift=0, middle=middle)
        b07otf = plot_with_rbs(b07otf_path, b07rbs_path, 'D', 5.0E-07, 0, r_shift=0, middle=middle)

        ignore = 0

        fig = pp.pplot(b07itf['LAMS Romp'][ignore:], b07itf['LAMS W'][ignore:], fmt='-',  color=6, lw=3)
        fig = pp.pplot(b07otf['LAMS Romp'][ignore:], b07otf['LAMS W'][ignore:], fmt='--', color=6, lw=3, fig=fig)
        fig = pp.pplot(b05itf['LAMS Romp'][ignore:], b05itf['LAMS W'][ignore:], fmt='-',  color=8, lw=3, fig=fig)
        fig = pp.pplot(b05otf['LAMS Romp'][ignore:], b05otf['LAMS W'][ignore:], fmt='--', color=8, lw=3, fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)')

    elif probe == 'c':

        middle = 1.25

        c05itf = plot_with_rbs(c05itf_path, c05rbs_path, 'D', 5.0E-07, 0, color=8, middle=middle)
        c05otf = plot_with_rbs(c05otf_path, c05rbs_path, 'U', 5.0E-07, 0, color=8, middle=middle)
        c07itf = plot_with_rbs(c07itf_path, c07rbs_path, 'U', 5.0E-07, 0, r_shift=0, middle=middle)
        c07otf = plot_with_rbs(c07otf_path, c07rbs_path, 'D', 5.0E-07, 0, r_shift=0, middle=middle)

        ignore = 0

        fig = pp.pplot(c07itf['LAMS Romp'][ignore:], c07itf['LAMS W'][ignore:], fmt='-',  color=6, lw=3)
        fig = pp.pplot(c07otf['LAMS Romp'][ignore:], c07otf['LAMS W'][ignore:], fmt='--', color=6, lw=3, fig=fig)
        fig = pp.pplot(c05itf['LAMS Romp'][ignore:], c05itf['LAMS W'][ignore:], fmt='-',  color=8, lw=3, fig=fig)
        fig = pp.pplot(c05otf['LAMS Romp'][ignore:], c05otf['LAMS W'][ignore:], fmt='--', color=8, lw=3, fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)')

    elif probe == 'a_rev':

        middle = 2.0

        a08itf = plot_with_rbs(a08itf_path, a08rbs_path, 'D', 0.5E-06, 0, color=8, middle=middle)
        a08otf = plot_with_rbs(a08otf_path, a08rbs_path, 'U', 0.5E-06, 0, color=8, middle=middle)

        ignore = 15

        fig = pp.pplot(a08itf['LAMS Romp'][ignore:], a08itf['LAMS W'][ignore:], fmt='-',  color=8, lw=3, label='ITF')
        fig = pp.pplot(a08otf['LAMS Romp'][ignore:], a08otf['LAMS W'][ignore:], fmt='--', color=8, lw=3, label='OTF', fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)')

    elif probe == 'a_for':

        middle = 2.0

        a15itf = plot_with_rbs(a15itf_path, a15rbs_path, 'U', 5.015E-07, 0, r_shift=0, middle=middle, avg=False)
        a15otf = plot_with_rbs(a15otf_path, a15rbs_path, 'D', 5.015E-07, 0, r_shift=0, middle=middle, avg=False)

        ignore = 70

        fig = pp.pplot(a15itf['LAMS Romp'][ignore:], a15itf['LAMS W'][ignore:], fmt='-',  color=6, lw=7, label='ITF')
        fig = pp.pplot(a15otf['LAMS Romp'][ignore+35:], a15otf['LAMS W'][ignore+35:], fmt='-', color=6, lw=3, label='OTF', fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)')
