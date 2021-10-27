import sys
sys.path.append("/Users/zamperini/github/oedge")
import oedge_plots


ncpath = "/Users/zamperini/Documents/d3d_work/184527/d3d-184527-bkg-003.nc"
omfit_path = "/Users/zamperini/Documents/d3d_work/184527/omfit_excel_file.csv"
op = oedge_plots.OedgePlots(ncpath)

# For v3.
#op.create_ts_from_omfit([184527], [2000,5000], omfit_path=[omfit_path],
#  smooth=True, filter=True, div_shift=-0.02, core_shift=0.0, smooth_window=21)

# For v4.
op.create_ts_from_omfit([184527], [2000,5000], omfit_path=[omfit_path],
  smooth=True, filter=True, div_shift=-0.00, core_shift=0.02, smooth_window=21)
