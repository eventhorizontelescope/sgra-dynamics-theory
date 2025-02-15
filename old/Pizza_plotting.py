'''  Pizza plot code. 

The first function, constraint_plot, is for Sgr A* Paper 5. The other functions are alternate versions.
It takes pass_table as its primary imput. Pass table is a 4d array. 2 (mad or sane) x 5 (spin) x n (inclination =5 values?) x n (rhigh = 4 values).
0 corresponds to a pass, 1 and -1 correspond to failing (although there's different options... I believe 1 and -1 are used for mutli=False keyword). 
When multi=False, 1 Failure corresponds to too high vaue, and -1 corresponds to too low. Multi=True ignores this and does something else. There are a few different failure states. 

Credit: Michi Baubock and others. I've been told the current pizza plot code is finnicky. Edit with caution.

'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib import *
import matplotlib
from matplotlib.gridspec import GridSpec
from matplotlib import patches

from definitions import *
#from read_files import *

# Set plotting style parameters
style_path = './Felix_4.mplstyle'
plt.rcParams.update(plt.rcParamsDefault)
plt.style.use(style_path)
#cmap = matplotlib.cm.get_cmap('viridis')
cmap = matplotlib.cm.get_cmap('RdYlGn_r')
matplotlib.rcParams['hatch.linewidth'] = 0.5
matplotlib.rcParams['font.family'] = 'serif'
#matplotlib.rcParams['mathtext.rm'] = 'Times'
#matplotlib.rcParams['mathtext.it'] = 'italic'
#matplotlib.rcParams['mathtext.default'] = 'stix'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
#plt.rcParams.update({"text.usetex":True})

#def constraint_plot(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['blue', 'green', 'red']):   #v0
#def constraint_plot(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['blue', 'green', 'red']):   #v1
#def constraint_plot(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['#F5793A', '#0F2080', '#85C0F9']):  #v2
def constraint_plot(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['#F5793A', 'darkgreen', '#85C0F9']):  #v2
    label_fontsize = 10

    #Make the colors more transparent for increasing number of failed models
    clow = col.to_rgba(color_tab[0])
    color_low = (clow[0], clow[1], clow[2], 1.0)
    color_low_t = (clow[0], clow[1], clow[2], 1.0)#0.5
#    color_low = (clow[0], clow[1], clow[2], 0.3)

    chigh =  col.to_rgba(color_tab[2])
    color_high = (chigh[0], chigh[1], chigh[2], 1.0)
    color_high_t = (chigh[0], chigh[1], chigh[2], 1.0) #0.5
#    color_high = (chigh[0], chigh[1], chigh[2], 0.8)

    color_good = col.to_rgba(color_tab[1])
    color_good_t = (color_good[0], color_good[1], color_good[2], 1.0)#0.5

    c1 = col.to_rgba(color_tab[1])
    color_fail1 = (c1[0], c1[1], c1[2], 0.3)

    c2 =  col.to_rgba(color_tab[1])
    color_fail2 = (c2[0], c2[1], c2[2], 0.2)

    c3 =  col.to_rgba(color_tab[1])
    color_fail3 = (c3[0], c3[1], c3[2], 0.1)

    color_ERROR = col.to_rgba('orange')

    nrows = len(flux_values)
    ncols = len(spin_values)

    fig, axs = plt.subplots(nrows, ncols, subplot_kw=dict(projection="polar"), figsize=(ncols*2, nrows*2),
                            gridspec_kw = dict(hspace=-0.02, wspace = 0.2))

    rlines = np.linspace(0., 1., len(rhigh_values) + 1, endpoint = True)
    tlines = np.array([0, 20, 40, 60, 80, 90])
    tlines = np.pi*tlines/180.

    th, r = np.meshgrid(tlines, rlines)
    thwidth = np.diff(th)
    rheight = np.diff(r, axis=0)

    th = th[:-1,:-1].flatten()
    r  = r[:-1,:-1].flatten()
    thwidth = thwidth[:-1].flatten()
    rheight = rheight[:,:-1].flatten()

    rlabelpos = (rlines[1:] + rlines[:-1]) / 2.
    tlabelpos = (tlines[1:] + tlines[:-1]) / 2.

    areas = []

    for m in flux_values:
        for c, s in enumerate(spin_values):
            axs[m,c].set_ylim([0,1])

            direction = (m % 2) * 2 - 1
            axs[m,c].set_theta_direction(direction)
            axs[m,c].set_thetamin(0)
            axs[m,c].set_thetamax(90)
            axs[m,c].set_theta_offset(direction * -np.pi / 2.)

            axs[m,c].set_xticks(tlines, fontsize = label_fontsize)
            axs[m,c].set_xticklabels('', fontsize = label_fontsize)
            axs[m,c].set_xticks(tlabelpos, minor=True, fontsize = label_fontsize)
            axs[m,c].set_xticklabels([r"{0:0d}$\degree$".format(i) for i in inc_values], minor=True, fontsize = label_fontsize)
            axs[m,c].set_yticks(rlines, fontsize = label_fontsize)
            axs[m,c].set_yticklabels('', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels(rhigh_values, minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            else:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels('', minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            axs[m,c].yaxis.tick_right()


            value = pass_table[m, c, :, :].transpose()
            colors = [['k' for c in value[0]] for r in value]
            for i in range(len(value)):
                for j in range(len(value[0])):
                    #print("(i,j) = " + str((i,j)))
                    if value[i,j] == 1:
                        if multi:
                            colors[i][j] = color_fail1
                            areas = areas + ['Fail 1']
                        else:
                            colors[i][j] = color_high
                            if i == 1 or j == 1 or j == 3:
                                colors[i][j] = color_high#_t
                            #areas = areas + ['Too high']
                            areas = areas + ['Fail All']
                    elif value[i,j] == 0:
                        colors[i][j] = color_good
                        if i == 1 or j == 1 or j == 3:
                            colors[i][j] = color_good#_t
                        #areas = areas + ['Pass']
                        areas = areas + ['Pass All']
                    elif value[i,j] == -1:
                        colors[i][j] = color_low
                        if i == 1 or j == 1 or j == 3:
                            colors[i][j] = color_low#_t
                        #areas = areas + ['Too low']
                        areas = areas + ['Pass Some']
                    elif value[i,j] == 2:
                        colors[i][j] = color_fail2
                        areas = areas + ['Fail 2']
                    elif value[i,j] >= 3:
                        colors[i][j] = color_fail3
                        areas = areas + ['Fail 3+']
                    elif value[i,j] == -2:
                        colors[i][j] = color_ERROR
                        areas = areas + ['Missing!']

            colors = np.array(colors)
            colors = colors.reshape((len(inc_values)*len(rhigh_values),4))

            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', color=colors)
            axs[m,c].grid(visible = True, which = 'major', axis = 'both', zorder = 50, linewidth = 0.7,
                          color = 'k')

            if m == 0:
                axs[m,c].set_title(f"$a_* = {s}$", fontsize = label_fontsize, pad = 8.0)
            if c == 0:
                axs[m,c].set_ylabel('MAD' if m == 1 else 'SANE', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', fontsize = label_fontsize, fontstyle = 'italic')
            else:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', va = 'top', fontsize = label_fontsize, fontstyle = 'italic')

    areas = np.unique(areas)
    #areas = np.array(['Pass All', 'Pass Some', 'Fail All'])
    patch_dict = {'Too high': patches.Patch(facecolor = color_high, edgecolor = 'k'),
                  'Pass': patches.Patch(facecolor = color_good, edgecolor = 'k'),
                  'Pass All': patches.Patch(facecolor = color_good, edgecolor = 'k'),
                  'Too low': patches.Patch(facecolor = color_low, edgecolor = 'k'),
                  'Fail 2': patches.Patch(facecolor = color_fail2, edgecolor = 'k'),
                  'Fail 3+': patches.Patch(facecolor = color_fail3, edgecolor = 'k'),
                  'Fail 1': patches.Patch(facecolor = color_fail1, edgecolor = 'k'),
                  'Missing!': patches.Patch(facecolor = color_ERROR, edgecolor = 'k'),
                  'Fail All': patches.Patch(facecolor = color_high, edgecolor = 'k'),
                  'Pass Some': patches.Patch(facecolor = color_low, edgecolor = 'k')}
    order_dict = {'Too high': 1, 'Fail 1': 1, 'Pass': 0, 'Too low': 2, 'Fail 2': 3, 'Fail 3+': 4, 'Missing!': 5, 'Pass All':0, 'Pass Some':1, 'Fail All':2}
    ordering = [order_dict[a] for a in areas]
    ind = np.argsort(ordering)
    areas = areas[ind]
    leg_patches = [patch_dict[a] for a in areas]
    fig.legend(leg_patches, areas, loc = 'lower center', ncol = len(areas), bbox_to_anchor = (0.5, 0.02), fontsize = label_fontsize)

    plt.figtext(0.1,0.49, r'$R_{\rm high}$', fontsize = label_fontsize)
    plt.suptitle(figtitle, fontsize = 18, y = 1.05)
    plt.savefig(filename, bbox_inches = 'tight')
    if show:
#         plt.tight_layout()
        plt.show()
    plt.close(fig)

    return


###
### Alternate versions:
###

def constraint_plot_h(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['blue', 'green', 'green', 'xkcd:red', 'xkcd:maroon']):

    #Make the colors more transparent for increasing number of failed models
    clow = col.to_rgba(color_tab[0])
    color_low = (clow[0], clow[1], clow[2], 0.3)
    chigh =  col.to_rgba(color_tab[1])
    color_high = (chigh[0], chigh[1], chigh[2], 0.3)
    color_good = col.to_rgba(color_tab[2])
    c2 =  col.to_rgba(color_tab[1])
    color_fail2 = (c2[0], c2[1], c2[2], 0.2)
    c3 =  col.to_rgba(color_tab[1])
    color_fail3 = (c3[0], c3[1], c3[2], 0.1)

    nrows = len(flux_values_h)
    ncols = len(spin_values_h)

    fig, axs = plt.subplots(nrows, ncols, subplot_kw=dict(projection="polar"), figsize=(ncols*2, nrows*2), gridspec_kw = dict(hspace=-0.02, wspace = 0.2))

    rlines = np.linspace(0., 1., len(rhigh_values_h) + 1, endpoint = True)
    tlines = np.array([0, 20, 70, 90])
    tlines = np.pi*tlines/180.

    th, r = np.meshgrid(tlines, rlines)
    thwidth = np.diff(th)
    rheight = np.diff(r, axis=0)

    th = th[:-1,:-1].flatten()
    r  = r[:-1,:-1].flatten()
    thwidth = thwidth[:-1].flatten()
    rheight = rheight[:,:-1].flatten()

    rlabelpos = (rlines[1:] + rlines[:-1]) / 2.
    tlabelpos = (tlines[1:] + tlines[:-1]) / 2.

    areas = []

    for m in flux_values_h:
        for c, s in enumerate(spin_values_h):
            axs[m,c].set_ylim([0,1])

            direction = (m % 2) * 2 - 1
            axs[m,c].set_theta_direction(direction)
            axs[m,c].set_thetamin(0)
            axs[m,c].set_thetamax(90)
            axs[m,c].set_theta_offset(direction * -np.pi / 2.)

            axs[m,c].set_xticks(tlines)
            axs[m,c].set_xticklabels('')
            axs[m,c].set_xticks(tlabelpos, minor=True)
            axs[m,c].set_xticklabels([r"{0:0d}$^\circ$".format(i) for i in inc_values_h], minor=True)
            axs[m,c].set_yticks(rlines)
            axs[m,c].set_yticklabels('')

            if m == 0:
                axs[m,c].set_yticks(rlabelpos, minor=True)
                axs[m,c].set_yticklabels(rhigh_values_h, minor=True)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            else:
                axs[m,c].set_yticks(rlabelpos, minor=True)
                axs[m,c].set_yticklabels('', minor=True)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            axs[m,c].yaxis.tick_right()


            value = pass_table[m, c, :, :].transpose()
#             value = np.zeros_like(value)
#             value[0] = np.ones(len(value[0]))
#             print(str(value))
            colors = [['k' for c in value[0]] for r in value]
            for i in range(len(value)):
                for j in range(len(value[0])):
                    if value[i,j] == 1:
                        if multi:
                            areas = areas + ['Fail 1']
                        else:
                            areas = areas + ['Too high']
                        colors[i][j] = color_high
                    elif value[i,j] == 0:
                        colors[i][j] = color_good
                        areas = areas + ['Pass']
                    elif value[i,j] == -1:
                        colors[i][j] = color_low
                        areas = areas + ['Too low']
                    elif value[i,j] == 2:
                        colors[i][j] = color_fail2
                        areas = areas + ['Fail 2']
                    elif value[i,j] >= 3:
                        colors[i][j] = color_fail3
                        areas = areas + ['Fail 3+']

            colors = np.array(colors)
            colors = colors.reshape((len(inc_values_h)*len(rhigh_values_h),4))

            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', color=colors)
            axs[m,c].grid(visible = True, which = 'major', axis = 'both', zorder = 50, linewidth = 0.7,
                          color = 'k')

            if m == 0:
                axs[m,c].set_title(f"$a = {s}$")
            if c == 0:
                axs[m,c].set_ylabel('MAD' if m == 1 else 'SANE', fontsize = 12)

            if m == 0:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$')
            else:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', va = 'top')

    areas = np.unique(areas)
    patch_dict = {'Too high': patches.Patch(facecolor = color_high),
                  'Pass': patches.Patch(facecolor = color_good),
                  'Too low': patches.Patch(facecolor = color_low),
                  'Fail 2': patches.Patch(facecolor = color_fail2),
                  'Fail 3+': patches.Patch(facecolor = color_fail3),
                  'Fail 1': patches.Patch(facecolor = color_high)}
    order_dict = {'Too high': 1, 'Fail 1': 1, 'Pass': 0, 'Too low': 2, 'Fail 2': 3, 'Fail 3+': 4}
    ordering = [order_dict[a] for a in areas]
    ind = np.argsort(ordering)
    areas = areas[ind]
    leg_patches = [patch_dict[a] for a in areas]
    fig.legend(leg_patches, areas, loc = 'lower center', ncol = len(areas), bbox_to_anchor = (0.5, 0.02))

    plt.figtext(0.1,0.49, r'$R_{\rm high}$')
    plt.suptitle(figtitle, fontsize = 18, y = 1.05)
    plt.savefig(filename, bbox_inches = 'tight')
    if show:
#         plt.tight_layout()
        plt.show()
    plt.close(fig)

    return

def constraint_plot_bhvk(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['blue', 'green', 'green', 'xkcd:red', 'xkcd:maroon']):

    #Make the colors more transparent for increasing number of failed models
    clow = col.to_rgba(color_tab[0])
    color_low = (clow[0], clow[1], clow[2], 0.3)
    chigh =  col.to_rgba(color_tab[1])
    color_high = (chigh[0], chigh[1], chigh[2], 0.3)
    color_good = col.to_rgba(color_tab[2])
    c2 =  col.to_rgba(color_tab[1])
    color_fail2 = (c2[0], c2[1], c2[2], 0.2)
    c3 =  col.to_rgba(color_tab[1])
    color_fail3 = (c3[0], c3[1], c3[2], 0.1)

    nrows = len(flux_values_bhvk)
    ncols = len(spin_values_bhvk)

    fig, axs = plt.subplots(nrows, ncols, subplot_kw=dict(projection="polar"), figsize=(ncols*2, nrows*2), gridspec_kw = dict(hspace=-0.02, wspace = 0.2))

    rlines = np.linspace(0., 1., len(rhigh_values_bhvk) + 1, endpoint = True)
    tlines = np.array([0, 20, 40, 60, 80, 90])
    tlines = np.pi*tlines/180.

    th, r = np.meshgrid(tlines, rlines)
    thwidth = np.diff(th)
    rheight = np.diff(r, axis=0)

    th = th[:-1,:-1].flatten()
    r  = r[:-1,:-1].flatten()
    thwidth = thwidth[:-1].flatten()
    rheight = rheight[:,:-1].flatten()

    rlabelpos = (rlines[1:] + rlines[:-1]) / 2.
    tlabelpos = (tlines[1:] + tlines[:-1]) / 2.

    areas = []

    for m in flux_values_bhvk:
        for c, s in enumerate(spin_values_bhvk):
            axs[m,c].set_ylim([0,1])

            direction = (m % 2) * 2 - 1
            axs[m,c].set_theta_direction(direction)
            axs[m,c].set_thetamin(0)
            axs[m,c].set_thetamax(90)
            axs[m,c].set_theta_offset(direction * -np.pi / 2.)

            axs[m,c].set_xticks(tlines)
            axs[m,c].set_xticklabels('')
            axs[m,c].set_xticks(tlabelpos, minor=True)
            axs[m,c].set_xticklabels([r"{0:0d}$^\circ$".format(i) for i in inc_values_bhvk], minor=True)
            axs[m,c].set_yticks(rlines)
            axs[m,c].set_yticklabels('')

            if m == 0:
                axs[m,c].set_yticks(rlabelpos, minor=True)
                axs[m,c].set_yticklabels(rhigh_values_bhvk, minor=True)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            else:
                axs[m,c].set_yticks(rlabelpos, minor=True)
                axs[m,c].set_yticklabels('', minor=True)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            axs[m,c].yaxis.tick_right()


            value = pass_table[m, c, :, :].transpose()
#             value = np.zeros_like(value)
#             value[0] = np.ones(len(value[0]))
#             print(str(value))
            colors = [['k' for c in value[0]] for r in value]
            for i in range(len(value)):
                for j in range(len(value[0])):
                    if value[i,j] == 1:
                        if multi:
                            areas = areas + ['Fail 1']
                        else:
                            areas = areas + ['Too high']
                        colors[i][j] = color_high
                    elif value[i,j] == 0:
                        colors[i][j] = color_good
                        areas = areas + ['Pass']
                    elif value[i,j] == -1:
                        colors[i][j] = color_low
                        areas = areas + ['Too low']
                    elif value[i,j] == 2:
                        colors[i][j] = color_fail2
                        areas = areas + ['Fail 2']
                    elif value[i,j] >= 3:
                        colors[i][j] = color_fail3
                        areas = areas + ['Fail 3+']

            colors = np.array(colors)
            colors = colors.reshape((len(inc_values_bhvk)*len(rhigh_values_bhvk),4))

            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', color=colors)
            axs[m,c].grid(visible = True, which = 'major', axis = 'both', zorder = 50, linewidth = 0.7,
                          color = 'k')

            if m == 0:
                axs[m,c].set_title(f"$a = {s}$")
            if c == 0:
                axs[m,c].set_ylabel('MAD' if m == 1 else 'SANE', fontsize = 12)

            if m == 0:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$')
            else:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', va = 'top')

    areas = np.unique(areas)
    patch_dict = {'Too high': patches.Patch(facecolor = color_high),
                  'Pass': patches.Patch(facecolor = color_good),
                  'Too low': patches.Patch(facecolor = color_low),
                  'Fail 2': patches.Patch(facecolor = color_fail2),
                  'Fail 3+': patches.Patch(facecolor = color_fail3),
                  'Fail 1': patches.Patch(facecolor = color_high)}
    order_dict = {'Too high': 1, 'Fail 1': 1, 'Pass': 0, 'Too low': 2, 'Fail 2': 3, 'Fail 3+': 4}
    ordering = [order_dict[a] for a in areas]
    ind = np.argsort(ordering)
    areas = areas[ind]
    leg_patches = [patch_dict[a] for a in areas]
    fig.legend(leg_patches, areas, loc = 'lower center', ncol = len(areas), bbox_to_anchor = (0.5, 0.02))

    plt.figtext(0.1,0.49, r'$R_{\rm high}$')
    plt.suptitle(figtitle, fontsize = 18, y = 1.05)
    plt.savefig(filename, bbox_inches = 'tight')
    if show:
#         plt.tight_layout()
        plt.show()
    plt.close(fig)

    return

def inverse_constraint_plot(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = True, color_tab = ['blue', 'green', 'green', 'xkcd:red', 'xkcd:maroon']):

    #Make the colors more transparent for increasing number of failed models
    c1 = col.to_rgba(color_tab[1])
    color1 = (c1[0], c1[1], c1[2], 1.0)
    c2 = col.to_rgba(color_tab[1])
    color2 = (c1[0], c1[1], c1[2], 0.7)
    c3 = col.to_rgba(color_tab[1])
    color3 = (c1[0], c1[1], c1[2], 0.4)
    c4 = col.to_rgba(color_tab[1])
    color4 = (c1[0], c1[1], c1[2], 0.1)
    c5 = col.to_rgba(color_tab[0])
    color5 = (c1[0], c1[1], c1[2], 0.0)

    nrows = len(flux_values)
    ncols = len(spin_values)

    fig, axs = plt.subplots(nrows, ncols, subplot_kw=dict(projection="polar"), figsize=(ncols*2, nrows*2),
                            gridspec_kw = dict(hspace=-0.02, wspace = 0.2))

    rlines = np.linspace(0., 1., len(rhigh_values) + 1, endpoint = True)
    tlines = np.array([0, 20, 40, 60, 80, 90])
    tlines = np.pi*tlines/180.

    th, r = np.meshgrid(tlines, rlines)
    thwidth = np.diff(th)
    rheight = np.diff(r, axis=0)

    th = th[:-1,:-1].flatten()
    r  = r[:-1,:-1].flatten()
    thwidth = thwidth[:-1].flatten()
    rheight = rheight[:,:-1].flatten()

    rlabelpos = (rlines[1:] + rlines[:-1]) / 2.
    tlabelpos = (tlines[1:] + tlines[:-1]) / 2.

    areas = []

    for m in flux_values:
        for c, s in enumerate(spin_values):
            axs[m,c].set_ylim([0,1])

            direction = (m % 2) * 2 - 1
            axs[m,c].set_theta_direction(direction)
            axs[m,c].set_thetamin(0)
            axs[m,c].set_thetamax(90)
            axs[m,c].set_theta_offset(direction * -np.pi / 2.)

            axs[m,c].set_xticks(tlines)
            axs[m,c].set_xticklabels('')
            axs[m,c].set_xticks(tlabelpos, minor=True)
            axs[m,c].set_xticklabels([r"{0:0d}$^\circ$".format(i) for i in inc_values], minor=True)
            axs[m,c].set_yticks(rlines)
            axs[m,c].set_yticklabels('')

            if m == 0:
                axs[m,c].set_yticks(rlabelpos, minor=True)
                axs[m,c].set_yticklabels(rhigh_values, minor=True)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            else:
                axs[m,c].set_yticks(rlabelpos, minor=True)
                axs[m,c].set_yticklabels('', minor=True)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            axs[m,c].yaxis.tick_right()


            value = pass_table[m, c, :, :].transpose()
#             value = np.zeros_like(value)
#             value[0] = np.ones(len(value[0]))
#             print(str(value))
            colors = [['k' for c in value[0]] for r in value]
            for i in range(len(value)):
                for j in range(len(value[0])):
                    if value[i,j] == 6:
                        colors[i][j] = color1
                        areas = areas + ['Fail 6']
                    elif value[i,j] == 5:
                        colors[i][j] = color2
                        areas = areas + ['Fail 5']
                    elif value[i,j] == 4:
                        colors[i][j] = color3
                        areas = areas + ['Fail 4']
                    elif value[i,j] == 3:
                        colors[i][j] = color4
                        areas = areas + ['Fail 3']
                    elif value[i,j] <= 2:
                        colors[i][j] = color5
                        areas = areas + ['Fail <= 2']

            colors = np.array(colors)
            colors = colors.reshape((len(inc_values)*len(rhigh_values),4))

            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', color=colors)
            axs[m,c].grid(visible = True, which = 'major', axis = 'both', zorder = 50, linewidth = 0.7,
                          color = 'k')

            if m == 0:
                axs[m,c].set_title(f"$a = {s}$")
            if c == 0:
                axs[m,c].set_ylabel('MAD' if m == 1 else 'SANE', fontsize = 12)

            if m == 0:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$')
            else:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', va = 'top')

    areas = np.unique(areas)
    patch_dict = {'Fail 6': patches.Patch(facecolor = color1),
                  'Fail 5': patches.Patch(facecolor = color2),
                  'Fail 4': patches.Patch(facecolor = color3),
                  'Fail 3': patches.Patch(facecolor = color4),
                  'Fail <= 2': patches.Patch(facecolor = color5)}
    order_dict = {'Fail 6': 0, 'Fail 5': 1, 'Fail 4': 2, 'Fail 3': 3, 'Fail <= 2': 4}
    ordering = [order_dict[a] for a in areas]
    ind = np.argsort(ordering)
    areas = areas[ind]
    leg_patches = [patch_dict[a] for a in areas]
    fig.legend(leg_patches, areas, loc = 'lower center', ncol = len(areas), bbox_to_anchor = (0.5, 0.02))

    plt.figtext(0.1,0.49, r'$R_{\rm high}$')
    plt.suptitle(figtitle, fontsize = 18, y = 1.05)
    plt.savefig(filename, bbox_inches = 'tight')
    if show:
#         plt.tight_layout()
        plt.show()
    plt.close(fig)

    return

def constraint_plot_paper_1(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['#F5793A', 'darkgreen', '#85C0F9']):  #v2

    label_fontsize = 10

    #Make the colors more transparent for increasing number of failed models
    clow = col.to_rgba(color_tab[0])
    color_low = (clow[0], clow[1], clow[2], 1.0)
    color_low_t = (clow[0], clow[1], clow[2], 1.0)#0.5
#    color_low = (clow[0], clow[1], clow[2], 0.3)

    chigh =  col.to_rgba(color_tab[2])
    color_high = (chigh[0], chigh[1], chigh[2], 1.0)
    color_high_t = (chigh[0], chigh[1], chigh[2], 1.0) #0.5
#    color_high = (chigh[0], chigh[1], chigh[2], 0.8)

    color_good = col.to_rgba(color_tab[1])
    color_good_t = (color_good[0], color_good[1], color_good[2], 1.0)#0.5

    c1 = col.to_rgba(color_tab[1])
    color_fail1 = (c1[0], c1[1], c1[2], 0.3)

    c2 =  col.to_rgba(color_tab[1])
    color_fail2 = (c2[0], c2[1], c2[2], 0.2)

    c3 =  col.to_rgba(color_tab[1])
    color_fail3 = (c3[0], c3[1], c3[2], 0.1)

    color_ERROR = col.to_rgba('orange')

    green = color_good
    orange = color_low
    red = color_high

    nrows = len(flux_values)
    ncols = len(spin_values)

    fig, axs = plt.subplots(nrows, ncols, subplot_kw=dict(projection="polar"), figsize=(ncols*2, nrows*2),
                            gridspec_kw = dict(hspace=-0.02, wspace = 0.2))

    rlines = np.linspace(0., 1., len(rhigh_values) + 1, endpoint = True)
    tlines = np.array([0, 20, 40, 60, 80, 90])
    tlines = np.pi*tlines/180.

    th, r = np.meshgrid(tlines, rlines)
    thwidth = np.diff(th)
    rheight = np.diff(r, axis=0)

    th = th[:-1,:-1].flatten()
    r  = r[:-1,:-1].flatten()
    thwidth = thwidth[:-1].flatten()
    rheight = rheight[:,:-1].flatten()

    rlabelpos = (rlines[1:] + rlines[:-1]) / 2.
    tlabelpos = (tlines[1:] + tlines[:-1]) / 2.

    areas = []

    for m in flux_values:
        for c, s in enumerate(spin_values):
            axs[m,c].set_ylim([0,1])

            direction = (m % 2) * 2 - 1
            axs[m,c].set_theta_direction(direction)
            axs[m,c].set_thetamin(0)
            axs[m,c].set_thetamax(90)
            axs[m,c].set_theta_offset(direction * -np.pi / 2.)

            axs[m,c].set_xticks(tlines, fontsize = label_fontsize)
            axs[m,c].set_xticklabels('', fontsize = label_fontsize)
            axs[m,c].set_xticks(tlabelpos, minor=True, fontsize = label_fontsize)
            axs[m,c].set_xticklabels([r"{0:0d}$^\circ$".format(i) for i in inc_values], minor=True, fontsize = label_fontsize)
            axs[m,c].set_yticks(rlines, fontsize = label_fontsize)
            axs[m,c].set_yticklabels('', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels(rhigh_values, minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            else:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels('', minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            axs[m,c].yaxis.tick_right()


            value = pass_table[m, c, :, :].transpose()
            colors = [['k' for c in value[0]] for r in value]
            for i in range(len(value)):
                for j in range(len(value[0])):
                    #print("(i,j) = " + str((i,j)))
                    #print(str(value[i,j]))
#                    if value[i,j] == 0:
#                        cval = 0.0
#                    else:
#                        cval = value[i,j] + 2.0
                    if value[i,j] == 0:
                        cval = 0.0
                    elif value[i,j] == 1 or value[i,j] == 2:
                        cval = 4
                    elif value[i,j] == 3:
                        cval = 5
                    elif  value[i,j] == 4 or value[i,j] == 5:
                        cval = 7
                    elif value[i,j] == 6 or value[i,j] == 7:
                        cval = 9
                    else:
                        cval = 10
                    colors[i][j] = cmap(cval / 10.0)
#                    if value[i,j] == 0:
#                        areas = areas + ['Pass All']
#                    if value[i,j] == 1:
#                        areas = areas + ['Some/All']
#                    if value[i,j] == 2:
#                        areas = areas + ['All/Some']
#                    if value[i,j] == 3:
#                        areas = areas + ['Some/Some']
#                    if value[i,j] == 4:
#                        areas = areas + ['None/All']
#                    if value[i,j] == 5:
#                        areas = areas + ['All/None']
#                    if value[i,j] == 6:
#                        areas = areas + ['None/Some']
#                    if value[i,j] == 7:
#                        areas = areas + ['Some/None']
#                    if value[i,j] == 8:
#                        areas = areas + ['None/None']

                    if value[i,j] == 0:
                        areas = areas + ['Pass All']
                    if value[i,j] == 1 or value[i,j] == 2:
                        areas = areas + ['Some/All']
                    if value[i,j] == 3:
                        areas = areas + ['Some/Some']
                    if value[i,j] == 4 or value[i,j] == 5:
                        areas = areas + ['All/None']
                    if value[i,j] == 6 or value[i,j] == 7:
                        areas = areas + ['Some/None']
                    if value[i,j] == 8:
                        areas = areas + ['None/None']

            colors = np.array(colors)
            colors = colors.reshape((len(inc_values)*len(rhigh_values),4))

            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', color=colors)
            axs[m,c].grid(visible = True, which = 'major', axis = 'both', zorder = 50, linewidth = 0.7,
                          color = 'k')

            if m == 0:
                axs[m,c].set_title(f"$a = {s}$", fontsize = label_fontsize, pad = 8.0)
            if c == 0:
                axs[m,c].set_ylabel('MAD' if m == 1 else 'SANE', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', fontsize = label_fontsize)
            else:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', va = 'top', fontsize = label_fontsize)

    areas = np.unique(areas)
    #areas = np.array(['Pass All', 'Pass Some', 'Fail All'])
    patch_dict = {'Pass All': patches.Patch(facecolor = cmap(0./10.), edgecolor = 'k'),
                  'Some/All': patches.Patch(facecolor = cmap(3./10.), edgecolor = 'k'),
                  'All/Some': patches.Patch(facecolor = cmap(4./10.), edgecolor = 'k'),
                  'Some/Some': patches.Patch(facecolor = cmap(5./10.), edgecolor = 'k'),
                  'None/All': patches.Patch(facecolor = cmap(6./10.), edgecolor = 'k'),
                  'All/None': patches.Patch(facecolor = cmap(7./10.), edgecolor = 'k'),
                  'None/Some': patches.Patch(facecolor = cmap(8./10.), edgecolor = 'k'),
                  'Some/None': patches.Patch(facecolor = cmap(9./10.), edgecolor = 'k'),
                  'None/None': patches.Patch(facecolor = cmap(10./10.), edgecolor = 'k')}
    order_dict = {'Pass All':  0,
                  'Some/All':  1,
                  'All/Some':  2,
                  'Some/Some': 3,
                  'None/All':  4,
                  'All/None':  5,
                  'None/Some': 6,
                  'Some/None': 7,
                  'None/None': 8}
    ordering = [order_dict[a] for a in areas]
    ind = np.argsort(ordering)
    areas = areas[ind]
    leg_patches = [patch_dict[a] for a in areas]
    fig.legend(leg_patches, areas, loc = 'lower center', ncol = len(areas), bbox_to_anchor = (0.5, 0.02), fontsize = label_fontsize, handletextpad = 0.1, handlelength = 1.0, columnspacing = 1.0)

    plt.figtext(0.1,0.49, r'$R_{\rm high}$', fontsize = label_fontsize)
    plt.suptitle(figtitle, fontsize = 18, y = 1.05)
    plt.savefig(filename, bbox_inches = 'tight')
    if show:
#         plt.tight_layout()
        plt.show()
    plt.close(fig)

    return

def constraint_plot_paper_1_hatched(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['#F5793A', 'darkgreen', '#85C0F9']):  #v2

    label_fontsize = 10

    #Make the colors more transparent for increasing number of failed models
    clow = col.to_rgba(color_tab[0])
    color_low = (clow[0], clow[1], clow[2], 1.0)
    color_low_t = (clow[0], clow[1], clow[2], 1.0)#0.5
#    color_low = (clow[0], clow[1], clow[2], 0.3)

    chigh =  col.to_rgba(color_tab[2])
    color_high = (chigh[0], chigh[1], chigh[2], 1.0)
    color_high_t = (chigh[0], chigh[1], chigh[2], 1.0) #0.5
#    color_high = (chigh[0], chigh[1], chigh[2], 0.8)

    color_good = col.to_rgba(color_tab[1])
    color_good_t = (color_good[0], color_good[1], color_good[2], 1.0)#0.5

    c1 = col.to_rgba(color_tab[1])
    color_fail1 = (c1[0], c1[1], c1[2], 0.3)

    c2 =  col.to_rgba(color_tab[1])
    color_fail2 = (c2[0], c2[1], c2[2], 0.2)

    c3 =  col.to_rgba(color_tab[1])
    color_fail3 = (c3[0], c3[1], c3[2], 0.1)

    color_ERROR = col.to_rgba('orange')

    green = color_good
    orange = color_low
    red = color_high

    nrows = len(flux_values)
    ncols = len(spin_values)

    fig, axs = plt.subplots(nrows, ncols, subplot_kw=dict(projection="polar"), figsize=(ncols*2, nrows*2),
                            gridspec_kw = dict(hspace=-0.02, wspace = 0.2))

    rlines = np.linspace(0., 1., len(rhigh_values) + 1, endpoint = True)
    tlines = np.array([0, 20, 40, 60, 80, 90])
    tlines = np.pi*tlines/180.

    th, r = np.meshgrid(tlines, rlines)
    thwidth = np.diff(th)
    rheight = np.diff(r, axis=0)

    th = th[:-1,:-1].flatten()
    r  = r[:-1,:-1].flatten()
    thwidth = thwidth[:-1].flatten()
    rheight = rheight[:,:-1].flatten()

    rlabelpos = (rlines[1:] + rlines[:-1]) / 2.
    tlabelpos = (tlines[1:] + tlines[:-1]) / 2.

    areas = []

    for m in flux_values:
        for c, s in enumerate(spin_values):
            axs[m,c].set_ylim([0,1])

            direction = (m % 2) * 2 - 1
            axs[m,c].set_theta_direction(direction)
            axs[m,c].set_thetamin(0)
            axs[m,c].set_thetamax(90)
            axs[m,c].set_theta_offset(direction * -np.pi / 2.)

            axs[m,c].set_xticks(tlines, fontsize = label_fontsize)
            axs[m,c].set_xticklabels('', fontsize = label_fontsize)
            axs[m,c].set_xticks(tlabelpos, minor=True, fontsize = label_fontsize)
            axs[m,c].set_xticklabels([r"{0:0d}$^\circ$".format(i) for i in inc_values], minor=True, fontsize = label_fontsize)
            axs[m,c].set_yticks(rlines, fontsize = label_fontsize)
            axs[m,c].set_yticklabels('', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels(rhigh_values, minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            else:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels('', minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            axs[m,c].yaxis.tick_right()


            value = pass_table[m, c, :, :].transpose()
            colors = [[(1.,1.,1.,0.) for c in value[0]] for r in value]
            int_hatchcolors = [[(0., 0., 0., 1.) for c in value[0]] for r in value]
            int_hatches = [[None for c in value[0]] for r in value]
            non_int_hatchcolors = [[(0., 0., 0., 1.) for c in value[0]] for r in value]
            non_int_hatches = [[None for c in value[0]] for r in value]
            for i in range(len(value)):
                for j in range(len(value[0])):
                    #print("(i,j) = " + str((i,j)))
                    #print(str(value[i,j]))
                    if value[i,j] == 0:
                        colors[i][j] = green
                    if value[i,j] == 1:
                        int_hatches[i][j] = '//'
                        int_hatchcolors[i][j] = green
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = orange
                    if value[i,j] == 2:
                        int_hatches[i][j] = '//'
                        int_hatchcolors[i][j] = orange
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = green
                    if value[i,j] == 3:
                        colors[i][j] = orange
                    if value[i,j] == 4:
                        int_hatches[i][j] = '//'
                        int_hatchcolors[i][j] = green
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = red
                    if value[i,j] == 5:
                        int_hatches[i][j] = '//'
                        int_hatchcolors[i][j] = red
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = green
                    if value[i,j] == 6:
                        int_hatches[i][j] = '//'
                        int_hatchcolors[i][j] = orange
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = red
                    if value[i,j] == 7:
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = orange
                        int_hatches[i][j] = '//'
                        int_hatchcolors[i][j] = red
                    if value[i,j] == 8:
                        colors[i][j] = red


                    if value[i,j] == 0:
                        areas = areas + ['Pass All']
                    if value[i,j] == 1 or value[i,j] == 2 or value[i,j] == 3 or value[i,j] == 6 or value[i,j] == 7:
                        areas = areas + ['Pass Some']
                    if value[i,j] == 4 or value[i,j] == 5 or value[i,j] == 8:
                        areas = areas + ['Fail All']

            colors = np.array(colors)
            colors = colors.reshape((len(inc_values)*len(rhigh_values),4))

            int_hatchcolors = np.array(int_hatchcolors)
            int_hatchcolors = int_hatchcolors.reshape((len(inc_values)*len(rhigh_values),4))
            int_hatches = np.array(int_hatches)
            int_hatches = int_hatches.reshape((len(inc_values)*len(rhigh_values)))

            non_int_hatchcolors = np.array(non_int_hatchcolors)
            non_int_hatchcolors = non_int_hatchcolors.reshape((len(inc_values)*len(rhigh_values),4))
            non_int_hatches = np.array(non_int_hatches)
            non_int_hatches = non_int_hatches.reshape((len(inc_values)*len(rhigh_values)))

            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', edgecolor = int_hatchcolors, color = colors, hatch = int_hatches)
            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', edgecolor = non_int_hatchcolors, color = colors, hatch = non_int_hatches)
            #axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', edgecolor = 'k', color = (0., 0., 0., 0.))#color=colors)
            axs[m,c].grid(visible = True, which = 'major', axis = 'both', zorder = 50, linewidth = 0.7,
                          color = 'k')
            for i,p in enumerate(axs[m,c].patches):
                p.set

            if m == 0:
                axs[m,c].set_title(f"$a = {s}$", fontsize = label_fontsize, pad = 8.0)
            if c == 0:
                axs[m,c].set_ylabel('MAD' if m == 1 else 'SANE', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', fontsize = label_fontsize)
            else:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', va = 'top', fontsize = label_fontsize)

    areas = areas + ['EHT', 'non-EHT']
    areas = np.unique(areas)
    #areas = np.array(['Pass All', 'Pass Some', 'Fail All'])
    patch_dict = {'Pass All': patches.Patch(facecolor = green, edgecolor = 'k'),
                  'Pass Some': patches.Patch(facecolor = orange, edgecolor = 'k'),
                  'Fail All': patches.Patch(facecolor = red, edgecolor = 'k'),
                  'EHT': patches.Patch(facecolor = (1., 1., 1., 0.), edgecolor = 'k', hatch = '//'),
                  'non-EHT': patches.Patch(facecolor = (1., 1., 1., 0.), edgecolor = 'k', hatch = '\\\\')}
    order_dict = {'Pass All':  0,
                  'Pass Some':  1,
                  'Fail All':  2,
                  'EHT': 3,
                  'non-EHT': 4}
    ordering = [order_dict[a] for a in areas]
    ind = np.argsort(ordering)
    areas = areas[ind]
    leg_patches = [patch_dict[a] for a in areas]
    fig.legend(leg_patches, areas, loc = 'lower center', ncol = len(areas), bbox_to_anchor = (0.5, 0.02), fontsize = label_fontsize, handlelength = 4.0)

    plt.figtext(0.1,0.49, r'$R_{\rm high}$', fontsize = label_fontsize)
    plt.suptitle(figtitle, fontsize = 18, y = 1.05)
    plt.savefig(filename, bbox_inches = 'tight')
    if show:
#         plt.tight_layout()
        plt.show()
    plt.close(fig)

    return

def constraint_plot_paper_1_hatched_v2(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['#F5793A', 'darkgreen', '#85C0F9']):  #v2

    label_fontsize = 10

    #Make the colors more transparent for increasing number of failed models
    clow = col.to_rgba(color_tab[0])
    color_low = (clow[0], clow[1], clow[2], 1.0)
    color_low_t = (clow[0], clow[1], clow[2], 1.0)#0.5
#    color_low = (clow[0], clow[1], clow[2], 0.3)

    chigh =  col.to_rgba(color_tab[2])
    color_high = (chigh[0], chigh[1], chigh[2], 1.0)
    color_high_t = (chigh[0], chigh[1], chigh[2], 1.0) #0.5
#    color_high = (chigh[0], chigh[1], chigh[2], 0.8)

    color_good = col.to_rgba(color_tab[1])
    color_good_t = (color_good[0], color_good[1], color_good[2], 1.0)#0.5

    c1 = col.to_rgba(color_tab[1])
    color_fail1 = (c1[0], c1[1], c1[2], 0.3)

    c2 =  col.to_rgba(color_tab[1])
    color_fail2 = (c2[0], c2[1], c2[2], 0.2)

    c3 =  col.to_rgba(color_tab[1])
    color_fail3 = (c3[0], c3[1], c3[2], 0.1)

    color_ERROR = col.to_rgba('orange')

    green = color_good
    orange = color_low
    red = color_high

    nrows = len(flux_values)
    ncols = len(spin_values)

    fig, axs = plt.subplots(nrows, ncols, subplot_kw=dict(projection="polar"), figsize=(ncols*2, nrows*2),
                            gridspec_kw = dict(hspace=-0.02, wspace = 0.2))

    rlines = np.linspace(0., 1., len(rhigh_values) + 1, endpoint = True)
    tlines = np.array([0, 20, 40, 60, 80, 90])
    tlines = np.pi*tlines/180.

    th, r = np.meshgrid(tlines, rlines)
    thwidth = np.diff(th)
    rheight = np.diff(r, axis=0)

    th = th[:-1,:-1].flatten()
    r  = r[:-1,:-1].flatten()
    thwidth = thwidth[:-1].flatten()
    rheight = rheight[:,:-1].flatten()

    rlabelpos = (rlines[1:] + rlines[:-1]) / 2.
    tlabelpos = (tlines[1:] + tlines[:-1]) / 2.

    areas = []

    for m in flux_values:
        for c, s in enumerate(spin_values):
            axs[m,c].set_ylim([0,1])

            direction = (m % 2) * 2 - 1
            axs[m,c].set_theta_direction(direction)
            axs[m,c].set_thetamin(0)
            axs[m,c].set_thetamax(90)
            axs[m,c].set_theta_offset(direction * -np.pi / 2.)

            axs[m,c].set_xticks(tlines, fontsize = label_fontsize)
            axs[m,c].set_xticklabels('', fontsize = label_fontsize)
            axs[m,c].set_xticks(tlabelpos, minor=True, fontsize = label_fontsize)
            axs[m,c].set_xticklabels([r"{0:0d}$^\circ$".format(i) for i in inc_values], minor=True, fontsize = label_fontsize)
            axs[m,c].set_yticks(rlines, fontsize = label_fontsize)
            axs[m,c].set_yticklabels('', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels(rhigh_values, minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            else:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels('', minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            axs[m,c].yaxis.tick_right()


            value = pass_table[m, c, :, :].transpose()
            colors = [[(1.,1.,1.,0.) for c in value[0]] for r in value]
            int_hatchcolors = [[(0., 0., 0., 1.) for c in value[0]] for r in value]
            int_hatches = [[None for c in value[0]] for r in value]
            non_int_hatchcolors = [[(0., 0., 0., 1.) for c in value[0]] for r in value]
            non_int_hatches = [[None for c in value[0]] for r in value]
            for i in range(len(value)):
                for j in range(len(value[0])):
                    #print("(i,j) = " + str((i,j)))
                    #print(str(value[i,j]))
                    if value[i,j] == 0:
                        colors[i][j] = green
                    if value[i,j] == 1:
                        colors[i][j] = green
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = orange
                    if value[i,j] == 2:
                        colors[i][j] = orange
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = green
                    if value[i,j] == 3:
                        colors[i][j] = orange
                    if value[i,j] == 4:
                        colors[i][j] = green
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = red
                    if value[i,j] == 5:
                        colors[i][j] = red
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = green
                    if value[i,j] == 6:
                        colors[i][j] = orange
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = red
                    if value[i,j] == 7:
                        non_int_hatches[i][j] = '\\\\'
                        non_int_hatchcolors[i][j] = orange
                        colors[i][j] = red
                    if value[i,j] == 8:
                        colors[i][j] = red


                    if value[i,j] == 0:
                        areas = areas + ['Pass All']
                    if value[i,j] == 1 or value[i,j] == 2 or value[i,j] == 3 or value[i,j] == 6 or value[i,j] == 7:
                        areas = areas + ['Pass Some']
                    if value[i,j] == 4 or value[i,j] == 5 or value[i,j] == 8:
                        areas = areas + ['Fail All']

            colors = np.array(colors)
            colors = colors.reshape((len(inc_values)*len(rhigh_values),4))

            int_hatchcolors = np.array(int_hatchcolors)
            int_hatchcolors = int_hatchcolors.reshape((len(inc_values)*len(rhigh_values),4))
            int_hatches = np.array(int_hatches)
            int_hatches = int_hatches.reshape((len(inc_values)*len(rhigh_values)))

            non_int_hatchcolors = np.array(non_int_hatchcolors)
            non_int_hatchcolors = non_int_hatchcolors.reshape((len(inc_values)*len(rhigh_values),4))
            non_int_hatches = np.array(non_int_hatches)
            non_int_hatches = non_int_hatches.reshape((len(inc_values)*len(rhigh_values)))

            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', edgecolor = int_hatchcolors, color = colors)
            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', edgecolor = non_int_hatchcolors, color = colors, hatch = non_int_hatches)
            #axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', edgecolor = 'k', color = (0., 0., 0., 0.))#color=colors)
            axs[m,c].grid(visible = True, which = 'major', axis = 'both', zorder = 50, linewidth = 0.7,
                          color = 'k')
            for i,p in enumerate(axs[m,c].patches):
                p.set

            if m == 0:
                axs[m,c].set_title(f"$a = {s}$", fontsize = label_fontsize, pad = 8.0)
            if c == 0:
                axs[m,c].set_ylabel('MAD' if m == 1 else 'SANE', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', fontsize = label_fontsize)
            else:
                axs[m,c].text(45*np.pi/180., 1.2, r'$i$', va = 'top', fontsize = label_fontsize)

    areas = areas + ['EHT', 'non-EHT']
    areas = np.unique(areas)
    #areas = np.array(['Pass All', 'Pass Some', 'Fail All'])
    patch_dict = {'Pass All': patches.Patch(facecolor = green, edgecolor = 'k'),
                  'Pass Some': patches.Patch(facecolor = orange, edgecolor = 'k'),
                  'Fail All': patches.Patch(facecolor = red, edgecolor = 'k'),
                  'EHT': patches.Patch(facecolor = (0., 0., 0., 1.), edgecolor = 'k', hatch = None),
                  'non-EHT': patches.Patch(facecolor = (1., 1., 1., 0.), edgecolor = 'k', hatch = '\\\\')}
    order_dict = {'Pass All':  0,
                  'Pass Some':  1,
                  'Fail All':  2,
                  'EHT': 3,
                  'non-EHT': 4}
    ordering = [order_dict[a] for a in areas]
    ind = np.argsort(ordering)
    areas = areas[ind]
    leg_patches = [patch_dict[a] for a in areas]
    fig.legend(leg_patches, areas, loc = 'lower center', ncol = len(areas), bbox_to_anchor = (0.5, 0.02), fontsize = label_fontsize, handlelength = 4.0)

    plt.figtext(0.1,0.49, r'$R_{\rm high}$', fontsize = label_fontsize)
    plt.suptitle(figtitle, fontsize = 18, y = 1.05)
    plt.savefig(filename, bbox_inches = 'tight')
    if show:
#         plt.tight_layout()
        plt.show()
    plt.close(fig)

    return

def constraint_plot_paper_1_hatched_v3(pass_table, filename = 'Constraint_Plot.png', figtitle = '', show = False, multi = False, color_tab = ['#F5793A', 'darkgreen', '#85C0F9']):  #v2

    label_fontsize = 10

    #Make the colors more transparent for increasing number of failed models
    clow = col.to_rgba(color_tab[0])
    color_low = (clow[0], clow[1], clow[2], 1.0)
    color_low_t = (clow[0], clow[1], clow[2], 1.0)#0.5
#    color_low = (clow[0], clow[1], clow[2], 0.3)

    chigh =  col.to_rgba(color_tab[2])
    color_high = (chigh[0], chigh[1], chigh[2], 1.0)
    color_high_t = (chigh[0], chigh[1], chigh[2], 1.0) #0.5
#    color_high = (chigh[0], chigh[1], chigh[2], 0.8)

    color_good = col.to_rgba(color_tab[1])
    color_good_t = (color_good[0], color_good[1], color_good[2], 1.0)#0.5

    c1 = col.to_rgba(color_tab[1])
    color_fail1 = (c1[0], c1[1], c1[2], 0.3)

    c2 =  col.to_rgba(color_tab[1])
    color_fail2 = (c2[0], c2[1], c2[2], 0.2)

    c3 =  col.to_rgba(color_tab[1])
    color_fail3 = (c3[0], c3[1], c3[2], 0.1)

    color_ERROR = col.to_rgba('orange')

    green = color_good
    orange = color_low
    red = color_high

    nrows = len(flux_values)
    ncols = len(spin_values)

    fig, axs = plt.subplots(nrows, ncols, subplot_kw=dict(projection="polar"), figsize=(ncols*2, nrows*2),
                            gridspec_kw = dict(hspace=-0.02, wspace = 0.2))

    rlines = np.linspace(0., 1., len(rhigh_values) + 1, endpoint = True)
    tlines = np.array([0, 20, 40, 60, 80, 90])
    tlines = np.pi*tlines/180.

    th, r = np.meshgrid(tlines, rlines)
    thwidth = np.diff(th)
    rheight = np.diff(r, axis=0)

    th = th[:-1,:-1].flatten()
    r  = r[:-1,:-1].flatten()
    thwidth = thwidth[:-1].flatten()
    rheight = rheight[:,:-1].flatten()

    rlabelpos = (rlines[1:] + rlines[:-1]) / 2.
    tlabelpos = (tlines[1:] + tlines[:-1]) / 2.

    areas = []

    for m in flux_values:
        for c, s in enumerate(spin_values):
            axs[m,c].set_ylim([0,1])

            direction = (m % 2) * 2 - 1
            axs[m,c].set_theta_direction(direction)
            axs[m,c].set_thetamin(0)
            axs[m,c].set_thetamax(90)
            axs[m,c].set_theta_offset(direction * -np.pi / 2.)

            axs[m,c].set_xticks(tlines, fontsize = label_fontsize)
            axs[m,c].set_xticklabels('', fontsize = label_fontsize)
            axs[m,c].set_xticks(tlabelpos, minor=True, fontsize = label_fontsize)
            #axs[m,c].set_xticklabels([r"{0:0d}$^\circ$".format(i) for i in inc_values], minor=True, fontsize = label_fontsize)
            axs[m,c].set_xticklabels([r"{0:0d}$\degree$".format(i) for i in inc_values], minor=True, fontsize = label_fontsize)
            axs[m,c].set_yticks(rlines, fontsize = label_fontsize)
            axs[m,c].set_yticklabels('', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels(rhigh_values, minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            else:
                axs[m,c].set_yticks(rlabelpos, minor=True, fontsize = label_fontsize)
                axs[m,c].set_yticklabels('', minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', grid_alpha = 1)
            axs[m,c].yaxis.tick_right()


            value = pass_table[m, c, :, :].transpose()
            colors = [[(1.,1.,1.,0.) for c in value[0]] for r in value]
            int_hatchcolors = [[(0., 0., 0., 1.) for c in value[0]] for r in value]
            int_hatches = [[None for c in value[0]] for r in value]
            non_int_hatchcolors = [[(0.3, 0.3, 0.3, 1.0) for c in value[0]] for r in value]
            non_int_hatches = [[None for c in value[0]] for r in value]
            for i in range(len(value)):
                for j in range(len(value[0])):
                    #print("(i,j) = " + str((i,j)))
                    #print(str(value[i,j]))
                    if value[i,j] == 0:
                        colors[i][j] = green
                    if value[i,j] == 1:
                        colors[i][j] = green
                        non_int_hatches[i][j] = '//////'
                    if value[i,j] == 2:
                        colors[i][j] = orange
                        non_int_hatches[i][j] = None
                    if value[i,j] == 3:
                        colors[i][j] = orange
                        non_int_hatches[i][j] = '//////'
                    if value[i,j] == 4:
                        colors[i][j] = green
                        non_int_hatches[i][j] = 'XXXXXX'
                    if value[i,j] == 5:
                        colors[i][j] = red
                        non_int_hatches[i][j] = None
                    if value[i,j] == 6:
                        colors[i][j] = orange
                        non_int_hatches[i][j] = 'XXXXXX'
                    if value[i,j] == 7:
                        non_int_hatches[i][j] = '//////'
                        non_int_hatchcolors[i][j] = (0.1, 0.1, 0.1, 1.0)
                        colors[i][j] = red
                    if value[i,j] == 8:
                        colors[i][j] = red
                        non_int_hatches[i][j] = 'XXXXXX'
                        non_int_hatchcolors[i][j] = (0.1, 0.1, 0.1, 1.0)


            colors = np.array(colors)
            colors = colors.reshape((len(inc_values)*len(rhigh_values),4))

            int_hatchcolors = np.array(int_hatchcolors)
            int_hatchcolors = int_hatchcolors.reshape((len(inc_values)*len(rhigh_values),4))
            int_hatches = np.array(int_hatches)
            int_hatches = int_hatches.reshape((len(inc_values)*len(rhigh_values)))

            non_int_hatchcolors = np.array(non_int_hatchcolors)
            non_int_hatchcolors = non_int_hatchcolors.reshape((len(inc_values)*len(rhigh_values),4))
            non_int_hatches = np.array(non_int_hatches)
            non_int_hatches = non_int_hatches.reshape((len(inc_values)*len(rhigh_values)))

            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', edgecolor = int_hatchcolors, color = colors)
            axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', edgecolor = non_int_hatchcolors, color = colors, hatch = non_int_hatches)
            #axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', edgecolor = 'k', color = (0., 0., 0., 0.))#color=colors)
            axs[m,c].grid(visible = True, which = 'major', axis = 'both', zorder = 50, linewidth = 0.7,
                          color = 'k')
            for i,p in enumerate(axs[m,c].patches):
                p.set

            if m == 0:
                astar_text = axs[m,c].set_title(f"$a_*$ = {s}", fontsize = label_fontsize,pad = 8.0, fontstyle = 'normal')
                #print(str(astar_text.get_window_extent()))
                #astar_text.set_clip_box(
            if c == 0:
                axs[m,c].set_ylabel('MAD' if m == 1 else 'SANE', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].text(45*np.pi/180., 1.2, r'i', fontsize = label_fontsize, fontstyle = 'italic')
            else:
                axs[m,c].text(45*np.pi/180., 1.2, r'i', va = 'top', fontsize = label_fontsize, fontstyle = 'italic')

    areas = ['EHT Pass All', 'EHT Pass Some', 'EHT Fail All', 'non-EHT Pass All', 'non-EHT Pass Some', 'non-EHT Fail All']
    areas = np.unique(areas)
    #areas = np.array(['Pass All', 'Pass Some', 'Fail All'])
    patch_dict = {'EHT Pass All': patches.Patch(facecolor = green, edgecolor = 'k'),
                  'EHT Pass Some': patches.Patch(facecolor = orange, edgecolor = 'k'),
                  'EHT Fail All': patches.Patch(facecolor = red, edgecolor = 'k'),
                  'non-EHT Pass All': patches.Patch(facecolor = (1., 1., 1., 0.), edgecolor = 'k', hatch = None),
                  'non-EHT Pass Some': patches.Patch(facecolor = (1., 1., 1., 0.), edgecolor = (0.3, 0.3, 0.3, 1.0), hatch = '//////'),
                  'non-EHT Fail All': patches.Patch(facecolor = (1., 1., 1., 0.), edgecolor = (0.3, 0.3, 0.3, 1.0), hatch = 'XXXXXX')}
    order_dict = {'EHT Pass All':  0,
                  'EHT Pass Some': 1,
                  'EHT Fail All':  2,
                  'non-EHT Pass All':  3,
                  'non-EHT Pass Some': 4,
                  'non-EHT Fail All':  5}
    ordering = [order_dict[a] for a in areas]
    ind = np.argsort(ordering)
    areas = areas[ind]
    leg_patches = [patch_dict[a] for a in areas]
    fig.legend(leg_patches, areas, loc = 'lower center', ncol = len(areas), bbox_to_anchor = (0.5, 0.02), fontsize = label_fontsize, handlelength = 1.0, handletextpad = 0.1)

    plt.figtext(0.1,0.49, r'$R_{\rm high}$', fontsize = label_fontsize)
    plt.suptitle(figtitle, fontsize = 18, y = 1.05)
    plt.savefig(filename, bbox_inches = 'tight')
    if show:
#         plt.tight_layout()
        plt.show()
    plt.close(fig)

    return

## For m87. Darker shading corresponds to higher likelihood. Feed the function a like_table, which requires specific shape/format. 
def likelihood_plot(like_table, incs = 0, magfield = 0, filename = 'Constraint_Plot', 
                    figtitle = '', show = False, cmap = blue_cmap, threshold = 0.0): 
    like_table_normalized = np.copy(like_table)
    like_inds = np.where(like_table_normalized < threshold) 
    like_table_normalized[like_inds] = threshold
    like_table_normalized -= np.min(like_table_normalized)
    like_table_normalized /= np.max(like_table_normalized)
#     like_table_normalized = like_table_normalized[np.where(like_table_normalized != -1)].\
#                             reshape((2, 5, 2, 4, 2, 2))
#     like_table_normalized /= np.max(like_table_normalized)
#     like_table_normalized += np.min(like_table_normalized[np.nonzero(like_table_normalized)])/1.1
#     like_table_normalized = np.log(like_table_normalized)
#     like_table_normalized -= np.min(like_table_normalized)
#     like_table_normalized /= np.max(like_table_normalized)
    
    cmissing = col.to_rgba('#F5793A')
    color_missing = (cmissing[0], cmissing[1], cmissing[2], 1.0)
    color_fail = (1.0, 1.0, 1.0, 1.0)
    
    label_fontsize = 10

    #Make the colors more transparent for increasing number of failed models
    notfound = "Not Found"

    nrows = len(all_ms)
    ncols = len(all_spins)
    
    
    like_temp = np.copy(like_table_normalized[:,:,0,:,0])
    like_comb = np.zeros_like(like_temp)
    
    
    for m,ms in enumerate(all_ms):
        for c, s in enumerate(all_spins):
            if incs == 0:
                value = like_table_normalized[m, c, 0, :, :, magfield].transpose()
            else:
                value = like_table_normalized[m, c, 1, :, :, magfield].transpose()
            for i in range(len(value)):
                for j in range(len(value[0])):
                    like_comb[m,c,j,i] = value[i,j]

    plot_rlows = all_rlows

    fig, axs = plt.subplots(nrows, ncols, subplot_kw=dict(projection="polar"), 
                            figsize=(ncols*2, nrows*2), gridspec_kw = dict(hspace=-0.02, wspace=0.2))

    rlines = np.linspace(0., 1., len(all_rhighs) + 1, endpoint = True)
    tlines = np.array([0, 45, 90])
    tlines = np.pi*tlines/180.

    th, r = np.meshgrid(tlines, rlines)
    thwidth = np.diff(th)
    rheight = np.diff(r, axis=0)

    th = th[:-1,:-1].flatten()
    r  = r[:-1,:-1].flatten()
    thwidth = thwidth[:-1].flatten()
    rheight = rheight[:,:-1].flatten()

    rlabelpos = (rlines[1:] + rlines[:-1]) / 2.
    tlabelpos = (tlines[1:] + tlines[:-1]) / 2.
    
    areas = []
#     inc_labels = []
    rl_labels = []
    rh_labels = [r"{0:0.0f}".format(rh) for rh in all_rhighs]
    for rl in plot_rlows:
        rl_labels += [r"{0:0.0f}".format(rl)]
            

    for m,ms in enumerate(all_ms):
        for c, s in enumerate(all_spins):
            axs[m,c].set_ylim([0,1])

            direction = (m % 2) * 2 - 1
            axs[m,c].set_theta_direction(direction)
            axs[m,c].set_thetamin(0)
            axs[m,c].set_thetamax(90)
            axs[m,c].set_theta_offset(direction * -np.pi / 2.)

            axs[m,c].set_xticks(tlines)#, fontsize = label_fontsize)
            axs[m,c].set_xticklabels('', fontsize = label_fontsize)
            axs[m,c].set_xticks(tlabelpos, minor=True)#, fontsize = label_fontsize)
            axs[m,c].set_xticklabels(rl_labels, minor=True, fontsize = label_fontsize)
            
            
            axs[m,c].set_yticks(rlines)#, fontsize = label_fontsize)
            axs[m,c].set_yticklabels('', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].set_yticks(rlabelpos, minor=True)#, fontsize = label_fontsize)
                axs[m,c].set_yticklabels(rh_labels, minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', 
                                               grid_alpha = 1)
            else:
                axs[m,c].set_yticks(rlabelpos, minor=True)#, fontsize = label_fontsize)
                axs[m,c].set_yticklabels('', minor=True, fontsize = label_fontsize)
                axs[m,c].yaxis.set_tick_params(which = 'major', length = 0)
                axs[m,c].yaxis.set_tick_params(which = 'minor', length = 0, grid_color = 'k', 
                                               grid_alpha = 1)
            axs[m,c].yaxis.tick_right()

            value_comb = like_comb[m, c, :, :]#.transpose()
            
            colors = [['k' for c in value_comb[0]] for r in value_comb]
            hatches = [['-' for c in value_comb[0]] for r in value_comb]
#             print(str(hatches))
        
            for i in range(len(value_comb)):
                for j in range(len(value_comb[0])):
                    colors[i][j] = cmap(value_comb[i,j])
                    if value_comb[i,j] == 0:
                        colors[i][j] = color_fail
                    if value_comb[i,j] == -1:
                        print("ms: " + str(ms) + ", spin: " + str(spin) + ": " + str((i,j)))
                    
#                     print(str(value_comb[i,j]))
                    hatches[i][j] = '-'

            colors = np.array(colors)
#             print(str(colors.shape))
#             colors = colors.reshape((len(plot_incs)*len(all_rhighs),4))
            colors = colors.reshape((len(plot_rlows)*len(all_rhighs),4))
#             print(str(colors.shape))
            hatches = np.array(hatches)
#             print(str(hatches.shape))
#             hatches = hatches.reshape((len(plot_incs)*len(all_rhighs)))
            hatches = hatches.reshape((len(plot_rlows)*len(all_rhighs)))

            temp = axs[m,c].bar(th, height=rheight, width=thwidth, bottom=r, align='edge', color=colors)
            axs[m,c].grid(visible = True, which = 'major', axis = 'both', zorder = 50, linewidth=1.2, 
                          color = 'k')
            
#             for i, bar in enumerate(temp.patches):
#                 barcol = bar.get_facecolor()
#                 if barcol == color_1:
#                     bar.set_edgecolor(color_2)
#                 if barcol == color_2:
#                     bar.set_edgecolor(color_1)

            if m == 0:
                axs[m,c].set_title(f"$a_* = {s}$", fontsize = label_fontsize, pad = 8.0)
            if c == 0:
                axs[m,c].set_ylabel('MAD' if m == 1 else 'SANE', fontsize = label_fontsize)

            if m == 0:
                axs[m,c].text(45*np.pi/180., 1.3, r'$R_{\rm low}$', fontsize = label_fontsize, 
                              fontstyle = 'italic')
            else:
                axs[m,c].text(45*np.pi/180., 1.3, r'$R_{\rm low}$', va = 'top', fontsize=label_fontsize,
                              fontstyle='italic')

    plt.figtext(0.1,0.49, r'$R_{\rm high}$', fontsize = label_fontsize)
    total_title = figtitle 
    plt.suptitle(total_title, fontsize = 18, y = 1.05)
    savename = filename + '_Hatched.png'
    plt.savefig(savename, bbox_inches = 'tight')
    if show:
        plt.show()
    plt.close(fig)

    return

