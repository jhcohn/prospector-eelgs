import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse


def add_sfh_plot(exout, fig, ax_loc=None,
                 main_color=None, tmin=False,  # tmin was 0.01
                 text_size=1, ax_inset=None, lw=1):
    '''
    add a small SFH plot at ax_loc
    text_size: multiply font size by this, to accomodate larger/smaller figures
    '''

    # set up plotting
    if ax_inset is None:
        if fig is None:
            ax_inset = ax_loc
        else:
            ax_inset = fig.add_axes(ax_loc, zorder=32)
    axfontsize = 4 * text_size

    xmin, ymin = np.inf, np.inf
    xmax, ymax = -np.inf, -np.inf
    for i, extra_output in enumerate(exout):

        #### load SFH
        t = extra_output['extras']['t_sfh']
        perc = np.zeros(shape=(len(t), 3))
        for jj in xrange(len(t)):
            perc[jj, :] = np.percentile(extra_output['extras']['sfh'][jj, :], [16.0, 50.0, 84.0])
            # 68.2% of population within 1 sigma <--> +/- 34.1%

        '''
        # REDSHIFT CONVERSION!
        z = []
        t = extra_output['extras']['t_sfh']
        for j in range(len(t)):
            z.append((2 / (3 * 7.22e-2 * t[j])) ** (2 / 3) - 1)
        '''
        #### plot SFH (t-->z)
        ax_inset.plot(t, perc[:, 1], '-', color=main_color[i], lw=lw)
        ax_inset.fill_between(t, perc[:, 0], perc[:, 2], color=main_color[i], alpha=0.3)
        ax_inset.plot(t, perc[:, 0], '-', color=main_color[i], alpha=0.3, lw=lw)
        ax_inset.plot(t, perc[:, 2], '-', color=main_color[i], alpha=0.3, lw=lw)

        #### update plot ranges
        # axis lims for nonparametric SFH
        xmin = np.min([xmin, t.min()])
        xmax = np.max([xmax, t.max()])
        ymin = np.min([ymin, perc.min()])
        ymax = np.max([ymax, perc.max()])

    #### labels, format, scales !
    if tmin:
        xmin = tmin

    # axlim_sfh = [xmax, xmin, ymin * .7, ymax * 1.4]
    axlim_sfh = [13.6, 10**-3, 10**-1, 10**3]
    # axlim_sfh = [4, 1., 10**-1, 10**3]  # redshift
    ax_inset.axis(axlim_sfh)
    ax_inset.set_ylabel(r'SFR [M$_{\odot}$/yr]', fontsize=axfontsize * 3, labelpad=2 * text_size)
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


# objs = [kwargs['obj']]
# for obj in objs:
def plotter(input, specific=False):

    # MAKE SURE OBJ AND FIELD ARE DEFINED
    obj = ''
    count = 0
    field = ''
    for i in input:
        if i == '_':
            count += 1
        elif count == 0:
            obj += i
        elif count == 1:
            field += i
        elif count == 2:
            break

    with open(input, 'rb') as file:
        extra_output = pickle.load(file)
        if specific:  # TESTING (print ssfr instead of regular sfr)
            plt.plot(extra_output['extras']['t_sfh'], extra_output['extras']['ssfr'], lw=2)
            plt.ylabel(r'Best-fit sSFR [yr$^{-1}$]')
        else:  # ORIGINAL
            plt.plot(extra_output['extras']['t_sfh'], extra_output['bfit']['sfh'], lw=2)
            plt.ylabel(r'Best-fit SFH [M$_\odot$ yr$^{-1}$]')
        plt.xlabel('t [Gyr]')
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rcParams.update({'font.size': 22})
        # plt.loglog()
        # ax = plt.gca()
        # ax.invert_xaxis()
        plt.title(field + '-' + obj)
        plt.show()

        if not specific:
            fig = plt.figure()
            sfh_ax = fig.add_axes([0.15, 0.15, 0.6, 0.6], zorder=32)
            add_sfh_plot([extra_output], fig, main_color=['black'], ax_inset=sfh_ax, text_size=3, lw=3)
            sfh_ax.invert_xaxis()
            plt.title(field + '-' + obj)
            plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--obj')
    parser.add_argument('--field')
    parser.add_argument('--spec')  # TESTING (including ssfr if --spec=True)

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    obj = kwargs['obj']
    field = kwargs['field']
    specific = kwargs['spec']

    if specific:
        base = '_sfh_out2.pkl'  # TESTING (includes ssfr)
    else:
        base = '_sfh_out.pkl'  # ORIGINAL
    file = obj + '_' + field + base

    plotter(file, specific=specific)  # TESTING: specific=True; ORIGINAL: specific=False (or remove "specific" keyword)

'''
# run from command line using (in snow environment): python print_sfh.py --obj=1824

python print_sfh.py --obj=1824 --field=cosmos

python print_sfh.py --obj=1824 --field=cosmos --spec=True  # TESTING
'''
