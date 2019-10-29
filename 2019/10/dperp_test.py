import oedge_plots
import pretty_plots as pp


nc003 = '/home/shawn/Documents/d3d_work/DIVIMP Runs/archive/Z0-167196-003.nc'
nc075 = '/home/shawn/Documents/d3d_work/DIVIMP Runs/Z0-167196-075.nc'
nc076 = '/home/shawn/Documents/d3d_work/DIVIMP Runs/Z0-167196-076b.nc'
nc077 = '/home/shawn/Documents/d3d_work/DIVIMP Runs/Z0-167196-077.nc'

op003 = oedge_plots.OedgePlots(nc003)  # 0.3
op075 = oedge_plots.OedgePlots(nc075)  # 1.0
op076 = oedge_plots.OedgePlots(nc076)  # 10.0
op077 = oedge_plots.OedgePlots(nc077)  # 0.03

x003, y003 = op003.along_ring(22, 'DDLIMS', charge='all', plot_it=False)
x075, y075 = op075.along_ring(22, 'DDLIMS', charge='all', plot_it=False)
x076, y076 = op076.along_ring(22, 'DDLIMS', charge='all', plot_it=False)
x077, y077 = op077.along_ring(22, 'DDLIMS', charge='all', plot_it=False)

fig = pp.pplot(x077, y077, label='0.03 m2/s', color=6,  fmt='-')
fig = pp.pplot(x003, y003, label='0.3  m2/s', color=8,  fmt='-', fig=fig)
fig = pp.pplot(x075, y075, label='1.0  m2/s', color=10, fmt='-', fig=fig)
fig = pp.pplot(x076, y076, label='10.0 m2/s', color=12, fmt='-', fig=fig,
               xlabel='s (m)', ylabel='W Density (m-3)', logy=True, xrange=[20, 70],
               yrange=[5e10, 1e16])
