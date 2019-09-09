import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


# These are the "Tableau 20" colors as RGB. The last one is just black. I added it.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
             (0, 0, 0)]

# A nice looking font.
plt.rcParams['font.family'] = 'serif'

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

def pplot(x, y, fmt='o', xerr=None, yerr=None, xlabel=None, ylabel=None, xrange=None,
          yrange=None, label=None, show_fig=True, alpha=1.0, weight='normal',
          fig=None, ms=8, fontsize=26, color=6, lw=5, logx=False, logy=False, linestyle=None):
    """
    A pretty plot for (x,y) data. Can be used to create just
    a single plot, or can return the figure instance if the user wants to
    make some extra changes beyond the default here. This is made to be
    a good starting point to create consistency among plots.

    To create multiple plots on one figure, call this function then pass the
    returned figure back into it with the other set(s) of data one at a time.
    The final figure should be what you want.
    """

#    plt.rcParams['font.family'] = 'serif'

    # These are the "Tableau 20" colors as RGB.
#    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
#                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
#                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
#                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
#                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
#    for i in range(len(tableau20)):
#        r, g, b = tableau20[i]
#        tableau20[i] = (r / 255., g / 255., b / 255.)

    # A good sized figure (4:3). Create a new one if not using iteratively.
    if fig is None:
        fig = plt.figure(figsize=(10, 7.5))
        ax1 = fig.add_subplot(111)
    else:
        ax1 = fig.axes[0]

    # Remove frame lines.
    ax1.spines["top"].set_visible(False)
    #ax1.spines["bottom"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    #ax1.spines["left"].set_visible(False)
    ax1.set_facecolor('white')

    # Axis ticks only on bottom and left.
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    # Set limits, if entered.
    if xrange is not None:
        try:
            ax1.set_xlim(xrange)
        except:
            print("Error: Make sure xrange is in the form [xmin, xmax].")
    if yrange is not None:
        try:
            ax1.set_ylim(yrange)
        except:
            print("Error: Make sure yrange is in the form [ymin, ymax].")

    # Make sure ticks are large enough to read.
    ax1.tick_params(axis='both', which='both', labelsize=18)
    ax1.set_xlabel(xlabel, fontsize=fontsize, weight=weight)
    ax1.set_ylabel(ylabel, fontsize=fontsize, weight=weight)

    # Set logarithmic scales if desired.
    if logx:
        plt.xscale('log')
    if logy:
        plt.yscale('log')
        #ax1.semilogx(x, y, fmt, label=label, ms=ms, color=tableau20[color], markeredgecolor='k', markeredgewidth='1', lw=lw)

    # Plot the data. Errorbar is more general in this case, so we use it.
    ax1.errorbar(x, y, fmt=fmt, xerr=xerr, yerr=yerr, label=label, ms=ms,
                 color=tableau20[color], markeredgecolor='k',
                 markeredgewidth='1', lw=lw, capsize=5, ecolor='k',
                 elinewidth=lw/5, capthick=lw/5, alpha=alpha, linestyle=linestyle)
    """
    ax1.plot(x, y, fmt, label=label, ms=ms, color=tableau20[color],
             markeredgecolor='k', markeredgewidth='1', lw=lw)
    """

    # Turn on the legend if a label was passed
    if label:
        ax1.legend(prop=dict(weight=weight, size=fontsize))

    if show_fig:
        fig.tight_layout()
        fig.show()

    return fig


def pcontourf(x, y, z, xlabel=None, ylabel=None, fontsize=26, weight='normal',
                cmap='plasma', cbarlabel='', extend='neither', xrange=None,
                yrange=None, vmin=None, vmax=None, show_fig=True):
    """
    Todo.
    """

    # A good sized figure.
    fig = plt.figure(figsize=(10, 7.5))
    ax1 = fig.add_subplot(111)

    # Remove frame lines.
    ax1.spines["top"].set_visible(False)
    #ax1.spines["bottom"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    #ax1.spines["left"].set_visible(False)
    ax1.set_facecolor('white')

    # Axis ticks only on bottom and left.
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    # Make sure ticks are large enough to read.
    ax1.tick_params(axis='both', which='both', labelsize=18)
    ax1.set_xlabel(xlabel, fontsize=fontsize, weight=weight)
    ax1.set_ylabel(ylabel, fontsize=fontsize, weight=weight)

    # Set limits, if entered.
    if xrange is not None:
        try:
            ax1.set_xlim(xrange)
        except:
            print("Error: Make sure xrange is in the form [xmin, xmax].")
    if yrange is not None:
        try:
            ax1.set_ylim(yrange)
        except:
            print("Error: Make sure yrange is in the form [ymin, ymax].")

    # Create countour plot and colorbar.
    cont = ax1.contourf(x, y, z, extend=extend, cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(cont)
    cbar.ax.set_ylabel(cbarlabel, size=fontsize)

    # Show the figure if desired.
    if show_fig:
        fig.tight_layout()
        fig.show()

    return fig

def pbar(y, y2=None, xlabel=None, ylabel=None, fontsize=26, color=6, color2=8, label=None,
         label2=None, weight='normal', bar_names=None, bar_names2=None, tick_rot=0, show_fig=True):

    # A good sized figure.
    fig = plt.figure(figsize=(10, 7.5))
    ax1 = fig.add_subplot(111)

    # Remove frame lines.
    ax1.spines["top"].set_visible(False)
    #ax1.spines["bottom"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    #ax1.spines["left"].set_visible(False)
    ax1.set_facecolor('white')

    # Axis ticks only on bottom and left.
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    # Make sure ticks are large enough to read.
    ax1.tick_params(axis='both', which='both', labelsize=18)
    ax1.set_xlabel(xlabel, fontsize=fontsize, weight=weight)
    ax1.set_ylabel(ylabel, fontsize=fontsize, weight=weight)

    # Create the bar graph.
    y_pos = np.arange(len(y))
    ax1.bar(y_pos, y, align='center', color=tableau20[color], label=label)

    # Add other bars if options added.
    if y2 is not None:
        y2_pos = np.arange(len(y2)) + len(y)
        bar_names = np.append(bar_names, bar_names2)
        ax1.bar(y2_pos, y2, align='center', color=tableau20[color2], label=label2)

    if bar_names is not None:
        ax1.set_xticks(np.arange(len(bar_names)))
        ax1.set_xticklabels(bar_names, rotation=tick_rot)

    # Turn on the legend if a label was passed
    if label:
        ax1.legend(prop=dict(weight=weight, size=fontsize))

    # Show the figure if desired.
    if show_fig:
        fig.tight_layout()
        fig.show()

    return fig

def show_colors():

    for i in range(0, len(tableau20)):
        plt.plot(np.linspace(0,0.1,10) * i, color=tableau20[i], label="Color "+str(i))
    plt.legend()
    plt.show()
