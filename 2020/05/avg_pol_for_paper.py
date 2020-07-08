import lim_plots as lim
import matplotlib.pyplot as plt


plt.rcParams['font.family'] = 'DejaVu Sans'
# These are the "Tableau 20" colors as RGB. The last one is just black. I added it.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
             (0, 0, 0)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

# Load the LimPlots objects for a simple and complex SOL case.
#simple_path  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-037.nc'
#complex_path = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-036.nc'
simple_path  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-049e.nc'
complex_path = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-050e.nc'
simp = lim.LimPlots(simple_path)
comp = lim.LimPlots(complex_path)

# Get the average poloidal data.
simp_pol = simp.avg_pol_profiles()
comp_pol = comp.avg_pol_profiles()

# Pull out some of the data into shorter variable names.
pol_locs = simp_pol['pol_locs'] * 100
simp_itf = simp_pol['avg_pol_itf']
simp_otf = simp_pol['avg_pol_otf']
comp_itf = comp_pol['avg_pol_itf']
comp_otf = comp_pol['avg_pol_otf']

fig = plt.figure()
ax = fig.add_subplot(111)

# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax.set_ylabel('Deposition (normalized)\n', fontsize=16)

ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
#ax1.plot(pol_locs, simp_itf/max(simp_itf.max(), comp_itf.max()), 'k', label='Simple')
#ax1.plot(pol_locs, comp_itf/max(simp_itf.max(), comp_itf.max()), 'k--', label='Complex')
#ax2.plot(pol_locs, simp_otf/max(simp_otf.max(), comp_otf.max()), 'k', label='Simple')
#ax2.plot(pol_locs, comp_otf/max(simp_otf.max(), comp_otf.max()), 'k--', label='Complex')
ax1.plot(pol_locs, simp_itf/simp_itf.max(), 'k', label='Simple')
ax1.plot(pol_locs, comp_itf/comp_itf.max(), 'k--', label='Complex')
ax2.plot(pol_locs, simp_otf/simp_otf.max(), 'k', label='Simple')
ax2.plot(pol_locs, comp_otf/comp_otf.max(), 'k--', label='Complex')
ax1.set_xlim([-0.25, 0.25])
ax2.set_xlim([-0.25, 0.25])
ax1.set_ylim([0.25, None])
ax2.set_ylim([0.25, None])
ax2.set_xlabel('Poloidal (cm)', fontsize=16)
ax2.legend(fontsize=16, loc='lower center')
#ax1.set_ylabel('Deposition', fontsize=16)
#ax2.set_ylabel('Deposition', fontsize=16)
fig.tight_layout()
fig.show()
