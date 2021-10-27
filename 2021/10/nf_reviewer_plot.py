import oedge_plots


ncpath = "/Users/zamperini/Documents/d3d_work/d3d-167247-bkg-013.nc"
op = oedge_plots.OedgePlots(ncpath)

pinato = op.read_data_2d("PINATO")
ne = op.read_data_2d("KNBS")

op.plot_contour_polygon("KTIBS", own_data=pinato/ne,
  cbar_label="neutral D / ne", vmin=0.0, vmax=0.25, lut=5)
