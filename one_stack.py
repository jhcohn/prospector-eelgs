import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import random
import os
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import uvj
'''
# preparing for printing priors  # MAYBE don't need; see if name=main
import eelg_fixedmet_params as ef_params
import eelg_thirty_params as et_params
import noelg_multirun_params as nm_params
import noelg_thirty_params as nt_params
'''


def draw_ssfr_from_prior(obj, fld, fig=None, axes=None, ndraw=1e4, alpha_sfh=1.0, pfile=None, show=True):
    # pfile=e_params or n_params

    # let's do it
    ndraw = int(ndraw)
    # zred = np.array([0.0, 0.5, 1.5, 2.5])  # where do we measure?
    zred = np.array([3.5])
    # want t_lookback in each bin I assume? Rather than these 4 zreds?
    logmass = np.array([10.])
    smass_factor = 0.8  # add in a not-unreasonable stellar mass <--> total mass conversion
    minssfr, maxssfr = 1e-14, 1e-5
    fs = 20

    # for i, z in enumerate(zred):
    # for index, redshift in zred: i.e. for (0, 0.0), (1, 0.5), (2, 1.5), (3, 2.5)

    # new redshift, new model
    # model = pfile.load_model(zred=z, alpha_sfh=alpha_sfh, **pfile.run_params)  # load_model(objname, field):
    model = pfile.load_model(objname=obj, field=fld)
    agebins = model.params['agebins']

    # create prior & draw from prior
    prior = model._config_dict['sfr_fraction']['prior']
    print(prior)
    mass = np.zeros(shape=(agebins.shape[0], ndraw))
    flatchain = np.random.dirichlet(tuple(1.0 for x in xrange(6)), ndraw)
    # flatchain = (array([0.stuff, ..., 0.stuff]), ndraw)  # array is ndraw arrays of 6 different 0.stuffs
    flatchain = flatchain[:, :-1]

    # 525 - 533
    # use fractions()  # need to convert logmass to linear mass before doing line 533
    for n in range(ndraw):
        mass[:, n] = pfile.sfrac_to_masses(logmass=logmass, agebins=agebins, sfrac=flatchain[n])

    time_per_bin = np.diff(10**agebins, axis=-1)[:, 0]

    ssfr = np.zeros(shape=(len(mass), len(mass[0])))  # 6, ndraw
    # print(len(ssfr[0]))
    # print(mass[0, :].sum(axis=0))
    for i in range(len(mass)):  # 6
        ssfr[i, :] = np.log10(mass[i, :] / time_per_bin[i] / 10**logmass)
        # ssfr[i, :] = np.log10(np.clip(mass[i, :].sum(axis=0) / time_per_bin[i].sum() / 10**logmass, minssfr, maxssfr))

    # print(ssfr)
    print('perc', np.percentile(ssfr[0], [16, 50, 84]))
    # median of each bin should have ssfr of 1/tuniv at that redshift (i.e. continuous SFH)
    # print(ssfr)  always [-7.]
    # print(len(mass[0]), len(mass))  # 10000, 6

    if show:
        # histogram
        fig = plt.figure()
        axes = plt.subplot(1, 1, 1)
        axes.hist(ssfr[-1], bins=50, histtype="step", normed=True)  # axes[i].hist()
        axes.set_xlim(np.log10(minssfr), np.log10(maxssfr))  # (-14, -7)  # axes[i].set_xlim()
        axes.set_ylim(0, 1)  # axes[i].set_ylim()

        axes.set_xlabel('sSFR (100 Myr)')
        # axes.set_xticks([])
        # axes.set_yticks([])
        '''
        if i > 1:
            axes[i].set_xlabel('sSFR (100 Myr)')
        else:
            axes[i].set_xticks([])

        if (i % 2) == 1:
            axes[i].set_yticks([])
        '''
        axes.text(0.02, 0.94, 'z = '+"{:.1f}".format(zred[0]), transform=axes.transAxes, fontsize=fs)  # .format(z)
        # axes[i].text(), axes[i].transAxes
        axes.text(0.02, 0.88, '<SFR>='+"{:.1f}".format(ssfr.mean()), transform=axes.transAxes, fontsize=fs)
        plt.show()
    else:
        return ssfr


def randraw(infile, num=1000):  # num=1000
    """
    For a given galaxy, randraw samples the posterior num times for each point in extra_output['extras']['sfh'][i]

    :param infile: ID_field_base_extra_out.py, where extra_output is stored using output.py
    :param num: number of times to sample the galaxy posterior at each point in the SFH
    :return: draw_from_sfh = 22 x num, lists the num random posterior samples, t_sfh = time vector associated with SFH
    """
    with open(infile, 'rb') as exout:
        extra_output = pickle.load(exout)

    draw_from_sfh = np.zeros(shape=(len(extra_output['extras']['ssfr']), num))  # shape=(22, num)
    # print(len(extra_output['extras']['ssfr']), len(extra_output['extras']['ssfr'][0]))  # 22, 2000
    # print(len(draw_from_sfh), len(draw_from_sfh[0]))  # 22, num

    for i in range(len(extra_output['extras']['ssfr'])):  # at each of these 22 points
        for j in range(num):  # randomly draw from the ssfr posterior num times
            draw_from_sfh[i][j] = extra_output['extras']['ssfr'][i][random.randint(0, num)]

    return draw_from_sfh, extra_output['extras']['t_sfh']


def smooth(perc):
    """
    Takes the stacked sfh that is output from stacker and averages the sfr values in each bin, such that the sfr within
    each bin is flat, as is the case in the original extra_output['extras']['sfh'] output

    :param perc: stored lists of the median and +/1 1sigma SFH values calculated from the gal_draws
    :return: smoother = same shape as perc, but with all points within a bin averaged s.t. all points within a bin=flat
    """
    # from perc: bin1 0:2, bin2 3:6, bin3 7:10, bin4 11:14, bin5 15:18, bin6 19:22
    smoother = np.zeros(shape=(len(perc), len(perc[0])))  # shape=(22, 3)
    for j in range(len(perc[0])):  # 3
        yax = perc[:, j]
        # print(yax)  # columns from perc: j=0 --> -1sigma, j=1 --> median, j=2 --> +1sigma
        for i in range(3):
            smoother[i][j] = (yax[0] + yax[1] + yax[2]) / 3
            smoother[i+19][j] = (yax[-1] + yax[-2] + yax[-3]) / 3
        for i in range(4):
            smoother[i+3][j] = (yax[3] + yax[4] + yax[5] + yax[6]) / 4
            smoother[i+7][j] = (yax[7] + yax[8] + yax[9] + yax[10]) / 4
            smoother[i+11][j] = (yax[11] + yax[12] + yax[13] + yax[14]) / 4
            smoother[i+15][j] = (yax[15] + yax[16] + yax[17] + yax[18]) / 4

    # print(smoother)
    return smoother


def stacker(gal_draws, sigma=1):
    """
    stacker takes input of random points drawn from a list of galaxies' SFH posteriors, concatenates them within each
    bin, and then calculates the median and 1 sigma errors in each bin
    gal_draws should be in format draws = [draw_from_sfh1, draw_from_sfh2, ...]
    each draw_from_sfh has shape=(22,num)

    :param gal_draws: list comprised of draw_from_sfh (i.e. comprised of the output from randraw) for a list of galaxies
    :param sigma: how many sigma of error we want to show in the plot
    :return: perc = stored lists of the median and +/- sigma SFH values calculated from the gal_draws
    """

    # len(gal_draws) = number of galaxies in stack; len(gal_draws[0]) = 22, len(gal_draws[0][0]) = num (1000)
    all_draws = np.zeros(shape=(len(gal_draws[0]), len(gal_draws[0][0]) * len(gal_draws)))
    for k in range(len(gal_draws)):
        # append the num=1000 values in each gal_draws[k] at each of the 22 points to all_draws:
        note = k * len(gal_draws[0][0])
        for i in range(len(gal_draws[k])):
            for j in range(len(gal_draws[k][i])):
                all_draws[i][note+j] += gal_draws[k][i][j]
    print(len(all_draws), len(all_draws[0]))  # 22, (number of galaxies in stack) * (num=1000)

    perc = np.zeros(shape=(len(gal_draws[0]), 2*sigma + 1))  # len(gal_draws[0])=22=len(t); len(perc)=22, len(perc[0])=3
    for jj in xrange(len(gal_draws[0])):
        if sigma == 1:
            perc[jj, :] = np.percentile(all_draws[jj, :], [16.0, 50.0, 84.0])  # median, +/- 34% = +/- 1sigma
        elif sigma == 3:
            perc[jj, :] = np.percentile(all_draws[jj, :], [0.3, 2.4, 16.0, 50.0, 84.0, 97.6, 99.7])  # out to 3 sigma

    return perc


def plot_sfhs(percs, t, lw=1, elist=None, uvj_in=False, spec=True, sigma=1, save=False, title=None,  priors=None,
              show=False):
    """
    Plots SFH stacks for two different galaxy samples side-by-side

    :param percs: list of two smoothed percs, for two different galaxy samples, each output by smooth(perc)
    :param t: time vector output by randraw
    :param lw: line width
    :param spec: if stacking specific SFR instead of plain SFR, spec=True
    :param sigma: how many sigma of error we want to show in the plot
    :return: plot
    """
    ymin, ymax = 1e-2, 1e3
    label = r'Stacked SFH [M$_\odot$ yr$^{-1}$]'

    if spec:
        ymin, ymax = 3e-11, 3e-8
        label = r'Stacked sSFH [yr$^{-1}$]'

    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)

    if uvj_in:  # also requires elist, llist to be not None; this insets uvj plots onto the top right of plot!
        inset_axes(ax1, width=3.2, height=2.8, loc=1)  # 20%
        uvj.uvj_plot(-1, 'all', objlist=elist, title=False, labels=False, lims=True, size=20, show=False)
        # create inset axis: width (%), height (inches), location
        # loc=1 (upper right), loc=2 (upper left) --> loc=3 (lower left), loc=4 (lower right); loc=7 (center right)
        # https://stackoverflow.com/questions/10824156/matplotlib-legend-location-numbers

    if sigma == 1:
        ax1.plot(t, percs[:, 1], '-', color='k', lw=lw)  # median
        ax1.fill_between(t, percs[:, 0], percs[:, 2], color='k', alpha=0.3)  # fill region between +/- 1sigma
        ax1.plot(t, percs[:, 0], '-', color='k', alpha=0.3, lw=lw)  # -1sigma
        ax1.plot(t, percs[:, 2], '-', color='k', alpha=0.3, lw=lw)  # +1sigma
        print('yo', percs[0, 1], percs[0, 2], percs[0, 0])

    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_ylim(ymin, ymax)  # , ymax
    ax1.set_xlim(10**-2, 2.5)  # (0, 2.5)  # (10**-2, 13.6)
    ax1.set_ylabel(label, fontsize=20)  # 30
    # ax1.text(4, 10**2.5, 'EELGs', fontsize=30)
    # ax1.text(1, 4*10**-8, 'EELGs', fontsize=30)  # use if not uvj_in

    if priors is not None:
        ssfr_pri = np.zeros(shape=(len(priors), 3))
        for i in range(len(priors)):  # 6
            ssfr_pri[i, :] = np.percentile(priors[i, :], [16, 50, 84])
        print(ssfr_pri)

        # new_t = [3.5*10**-2, 2.5*10**-1, 7*10**-1, 1.15, 1.4, 1.8]
        fill_t = [0, 100]
        y1 = [10 ** ssfr_pri[0][0]] * 2
        y2 = [10 ** ssfr_pri[0][2]] * 2
        ax1.axhline(y=10 ** ssfr_pri[0][1], color='r')
        ax1.axhline(y=10 ** ssfr_pri[0][0], color='r')
        ax1.axhline(y=10 ** ssfr_pri[0][2], color='r')
        ax1.fill_between(fill_t, y1, y2, color='r', hatch='/', facecolor='none')  # , alpha=0.3)

    plt.tight_layout()
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size': 20})
    ax1.set_xlabel(r'Lookback time [Gyr]', fontsize=20)
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=20)
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=20)
    # fig.text(0.5, 0.04, 'Lookback time [Gyr]', ha='center', fontsize=20)  # 30
    plt.tight_layout()
    if save:
        plt.savefig(title + '.png', bbox_inches='tight')
    else:
        plt.show()


if __name__ == "__main__":

    others = False
    if others:
        base = 'thirty'
        folders = 'etpkls/'
        import eelg_thirty_params as e_params
    else:
        base = 'fixedmet'  # use for fixedmet
        folders = 'pkls/'
        import eelg_fixedmet_params as e_params

    pri = draw_ssfr_from_prior('1824', 'cosmos', ndraw=1e4, alpha_sfh=1.0, pfile=e_params, show=False)

    eelg_list = open('eelg_specz_ids', 'r')
    pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders
    eelgs = []
    for line in eelg_list:
        if line[0] == '#':
            pass
        else:
            cols = line.split()
            eelgs.append(cols[1] + '_' + cols[0] + '_' + base)  # base[0] = fixedmet (or otherbins)
    eelg_list.close()

    # START STACKING
    t1 = []
    draws = []
    boots = []
    nummy = 0
    c = 0
    for glxy in eelgs:
        c += 1
        # file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
        file = pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            nummy += 1
            temp = randraw(file)  # temp[0] lists the num=1000 random posterior samples; temp[1] = time vector
            # temp = bootdraw(file)  # temp[0] lists the num=1000 random posterior samples; temp[1] = time vector
            draws.append(temp[0])
            # boots.append(bootstrap(temp[0]))
            t1.append(temp[1])
        else:
            print(file)

    # stacker2(draws, t1)
    sig = 1  # what sigma error to show on plot
    perc1 = stacker(draws, sigma=sig)

    print(nummy, c, 'nume')
    smooth_perc = smooth(perc1)
    plot_sfhs(smooth_perc, t1[0], elist=eelgs, uvj_in=True, sigma=sig, priors=pri)


'''
# run from command line in snow environment using:
python one_stack.py
'''
