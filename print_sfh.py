import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse


def add_sfh_plot(exout, fig, ax_loc=None, main_color=None, tmin=False, text_size=1,
                 ax_inset=None, lw=1, specific=False, thirty=True):
    """
    add a small SFH plot at ax_loc

    :param exout: extra_output generated and pickled by output.py
    :param fig:
    :param ax_loc:
    :param main_color:
    :param tmin:
    :param text_size: ultiply font size by this, to accommodate larger/smaller figures
    :param ax_inset:
    :param lw: line width
    :param specific: if True, plot sSFH
    :param thirty: if True, indicates I'm using different set of parameter files
    :return: SFH plot with 1sigma errors
    """

    # set up plotting
    if ax_inset is None:
        if fig is None:
            ax_inset = ax_loc
        else:
            ax_inset = fig.add_axes(ax_loc, zorder=32)
    axfontsize = 4 * text_size

    # set limits
    xmin, ymin = np.inf, np.inf
    xmax, ymax = -np.inf, -np.inf
    for i, extra_output in enumerate(exout):

        # load SFH from extra_output
        t = extra_output['extras']['t_sfh']

        perc = np.zeros(shape=(len(t), 3))  # initialize percentiles for SFH
        if specific:
            sfh = extra_output['extras']['ssfr']
            ymin, ymax = 10**-12, 10**-7
            ylab = r'SSFR [yr$^{-1}$]'
        else:
            sfh = extra_output['extras']['sfh']
            ymin, ymax = 10**-1, 10**3
            ylab = r'SFR [M$_{\odot}$/yr]'
        for jj in xrange(len(t)):
            perc[jj, :] = np.percentile(sfh[jj, :], [16.0, 50.0, 84.0])
            # print(extra_output['extras']['sfh'][jj, :], len(extra_output['extras']['sfh'][jj, :])) # ncalc (output.py)
            # 68.2% of population within 1 sigma <--> +/- 34.1% from the 50th percentile
            # np.percentile: array of requested percentiles for each point in SFH (or SSFH)
        # print(perc)

        # plot SFH
        ax_inset.plot(t, perc[:, 1], '-', color=main_color[i], lw=lw)  # plot median
        ax_inset.fill_between(t, perc[:, 0], perc[:, 2], color=main_color[i], alpha=0.3)  # fill region between +/-1sig
        ax_inset.plot(t, perc[:, 0], '-', color=main_color[i], alpha=0.3, lw=lw)  # plot -1sigma
        ax_inset.plot(t, perc[:, 2], '-', color=main_color[i], alpha=0.3, lw=lw)  # plot +1sigma

        # update plot ranges
        # axis lims for nonparametric SFH
        xmin = np.min([xmin, t.min()])
        xmax = np.max([xmax, t.max()])
        ymin = np.min([ymin, perc.min()])
        ymax = np.max([ymax, perc.max()])

    # labels, format, scales !
    if tmin:
        xmin = tmin
    elif thirty:
        xmin = 10 ** -3  # xmin at 10**-3 when using 30 Myr bin
    else:
        xmin = 10 ** -2  # xmin at 10**-2, unless using 10 Myr bin, in which case xmin at 10**-3

    # axlim_sfh = [13.6, xmin, 10**-1, 10**3]  # [xmax, xmin, ymin, ymax]
    axlim_sfh = [13.6, xmin, ymin, ymax]  # [xmax, xmin, ymin, ymax]
    # axlim_sfh = [4, 1., 10**-1, 10**3]  # redshift
    ax_inset.axis(axlim_sfh)
    ax_inset.set_ylabel(ylab, fontsize=axfontsize * 3, labelpad=2 * text_size)
    ax_inset.set_xlabel(r't$_{\mathrm{lookback}}$ [Gyr]', fontsize=axfontsize * 3, labelpad=2 * text_size)

    subsx = [2, 3, 4, 5, 6, 7, 8, 9]
    subsy = [2, 3, 4, 5, 6, 7, 8, 9]

    ax_inset.set_xscale('log', nonposx='clip', subsx=subsx)
    ax_inset.xaxis.set_tick_params(labelsize=axfontsize * 2)

    ax_inset.set_yscale('log', nonposy='clip', subsy=subsy)
    ax_inset.yaxis.set_tick_params(labelsize=axfontsize * 2)

    ax_inset.tick_params('both', length=lw * 3, width=lw * .6, which='major')

    ax_inset.tick_params('both', length=lw, width=lw * .3, which='minor')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax_inset.spines[axis].set_linewidth(lw * .6)


def plotter(input, specific=False, simple=False, thirty=False):
    """

    :param input: this is the pickled extra_output genereted by output.py
    :param specific: plot SFH if False, sSFH if True
    :param simple: if True, just plot best fit
    :param thirty: if True, indicates I'm using results from a different set of parameter files
    :return: SFH plot
    """
    # MAKE SURE OBJ AND FIELD ARE DEFINED
    obj = ''
    count = 0
    field = ''
    if input[0] == 'p':  # likely true for plotting in make_all_plots.py or all_seds.py
        slash = 0
        for i in input:
            if i == '/':
                slash += 1
            elif slash == 1:
                if i == '_':
                    count += 1
                elif count == 0:
                    obj += i
                elif count == 1:
                    field += i
                elif count == 2:
                    break
    else:
        for i in input:
            if i == '_':
                count += 1
            elif count == 0:
                obj += i
            elif count == 1:
                field += i
            elif count == 2:
                break

    plus = 0
    if field == 'cosmos':
        plus = 100000
    elif field == 'uds':
        plus = 200000
    print(int(obj) + plus)

    with open(input, 'rb') as file:
        extra_output = pickle.load(file)

        '''
        if specific:  # TESTING (print ssfr instead of regular sfr)
            plt.loglog(extra_output['extras']['t_sfh'], extra_output['extras']['ssfr'], lw=2, color='k')
            plt.ylabel(r'Best-fit sSFR [yr$^{-1}$]')
            plt.ylim(10**-13, 10**-7)
            plt.xlim(10**-2, 13.6)
        '''

        if simple:  # ORIGINAL
            plt.plot(extra_output['extras']['t_sfh'], extra_output['bfit']['sfh'], lw=2)
            plt.ylabel(r'Best-fit SFH [M$_\odot$ yr$^{-1}$]')
            plt.ylim(10**-2, 10**3)
            plt.xlim(10**-2, 13.6)

            plt.xlabel('t [Gyr]')
            plt.rc('xtick', labelsize=20)
            plt.rc('ytick', labelsize=20)

        else:  # see add_sfh_plot() function above!
            fig = plt.figure()
            sfh_ax = fig.add_axes([0.15, 0.15, 0.6, 0.6], zorder=32)
            add_sfh_plot([extra_output], fig, main_color=['black'], ax_inset=sfh_ax, text_size=3, lw=3,
                         specific=specific, thirty=thirty)
            sfh_ax.invert_xaxis()

        plt.rcParams.update({'font.size': 22})
        plt.title(field + '-' + obj)
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--obj')
    parser.add_argument('--field')
    parser.add_argument('--base')  # needs to have format --base=_newbins
    specific = False

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    obj = kwargs['obj']
    field = kwargs['field']
    base = kwargs['base']

    # select correct pkl based on file naming system
    file = obj + '_' + field + '_' + base + '_extra_out' + '.pkl'

    # plot SFH!
    plotter(file, specific=specific, simple=False)

'''
# REDSHIFT CONVERSION!
z = []
t = extra_output['extras']['t_sfh']
for j in range(len(t)):
    z.append((2 / (3 * 7.22e-2 * t[j])) ** (2 / 3) - 1)
'''
'''
# run from command line using:
python print_sfh.py --obj=1824 --field=cosmos --base=newbins
'''
