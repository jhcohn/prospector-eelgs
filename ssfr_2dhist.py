import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import random
import os
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import uvj
from astropy.cosmology import WMAP9
import scipy
from matplotlib.ticker import MultipleLocator, NullFormatter, MaxNLocator
from matplotlib import rc
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec


def draw_sfrs(objs, flds, new_logmass=None, ndraw=1e4, alpha_sfh=1.0, pfile=None, show=True):
    # pfile=e_params or n_params

    # let's do it
    ndraw = int(ndraw)
    zred = np.array([3.5])
    logmass = np.array([10.])
    smass_factor = 0.8  # add in a not-unreasonable stellar mass <--> total mass conversion
    minssfr, maxssfr = 1e-14, 1e-5

    # new redshift, new model
    if new_logmass is None:
        new_logmass = 10.17
    new_SFRs = []
    for i in range(len(objs)):
        new_model = pfile.load_model(objname=objs[i], field=flds[i])
        new_zred = new_model.params['zred']
        tuniv = WMAP9.age(new_zred).value
        new_SFRs.append(new_logmass / tuniv)

    return new_SFRs


def draw_ssfr_from_prior2(obj, fld, gals=None, ndraw=1e4, pfile=None):
    # let's do it
    ndraw = int(ndraw)
    if gals == 'eelgs':
        logmass = np.array([10.17])
    elif gals == 'lbgs':
        logmass = np.array([10.5])
    else:
        logmass = np.array([10.])
    # minssfr, maxssfr = 1e-14, 1e-5

    # new redshift, new model
    model = pfile.load_model(objname=obj, field=fld)  # load_model prints zred, tuniv
    agebins = model.params['agebins']

    # create prior & draw from prior
    prior = model._config_dict['sfr_fraction']['prior']
    print(prior, 'prior')

    mass = np.zeros(shape=(agebins.shape[0], ndraw))
    # create a prior and sample from it. what you want from this calculation is an array of masses in each bin drawn
    # directly from the prior. the two lines below will create a series of sfr_fractions drawn directly from the prior.
    # what you have to do is convert these sfr_fractions in 'flatchain' into masses in each bin, in an array with the
    # same shape as the 'mass' vector above. these calculations are done in lines 525-533 here:
    # https://github.com/jhcohn/prospector-eelgs/blob/master/eelg_fixedmet_params.py, where each line in 'flatchain'
    # below is a set of 'sfr_fractions'. Here you can assume the total log(M) = 10, as stated above.
    flatchain = np.random.dirichlet(tuple(1.0 for x in xrange(6)), ndraw)
    # flatchain = (array([0.stuff, ..., 0.stuff]), ndraw)
    flatchain = flatchain[:, :-1]  # array is ndraw arrays of 6 different 0.stuffs

    for n in range(ndraw):
        mass[:, n] = pfile.sfrac_to_masses(logmass=logmass, agebins=agebins, sfrac=flatchain[n])

    # convert to sSFR
    # print(mass.shape)  # (6, ndraw)
    # sprint(mass[0, :].sum(axis=0))  # consistently ~5.4e12
    # print(agebins)  # [[0, 8], [8, 8.7],...]
    time_per_bin = np.diff(10**agebins, axis=-1)[:, 0]
    # print(np.diff(10**agebins, axis=-1))  # [[10^8 - 10^0],[10^8.7 - 10^8],...]
    # print(time_per_bin)  # [1e8, 4.01e8, 4.99e8, 2.84e8, 3.65e8, 4.685e8]
    # print(time_per_bin[0].sum)  # "<built-in method sum of numpy.float64 object at 0x7f99da4993d8>"
    # print(mass[0, :].sum(axis=0) / time_per_bin[0].sum() / 10**logmass)  # 5.3e-6 to 5.4e-6 > 1e-7 --> clipped
    # print(np.clip(mass[0, :].sum(axis=0) / time_per_bin[0].sum() / 10**logmass, minssfr, maxssfr))  ~ 5e-6 or 1e-7

    # NOTE: time_per_bin[0] corresponds to 0-100 Myr time bin = most recent time bin
    ssfr = np.zeros(shape=(len(mass), len(mass[0])))  # 6, ndraw
    mass_cml = np.zeros(shape=(len(mass), len(mass[0])))

    for i in range(len(mass)):  # 6
        n = 0
        while n <= i:
            mass_cml[i, :] += mass[n, :]
            n += 1

    for i in range(len(mass)):
        # ssfr[i, :] = np.log10(mass[i, :] / time_per_bin[i] / mass_cml[-(i+1), :])  # mass_cml defined in reverse order
        # = (fractional mass formed in bin) / (time per bin) / (cumulative mass formed up to this point)
        ssfr[i, :] = np.log10(mass[i, :] / time_per_bin[i] / 10**logmass)  # mass_cml defined in reverse order
        # ssfr[i, :] = np.log10(np.clip(mass[i, :].sum(axis=0) / time_per_bin[i].sum() / 10**logmass, minssfr, maxssfr))

    perc = np.zeros(shape=(len(mass), 3))  # 6, 3
    for i in range(len(perc)):
        perc[i, :] = np.percentile(ssfr[i, :], [16., 50., 84.])
    # print('perc', np.percentile(ssfr[0], [16, 50, 84]))
    print(perc, 'perc')

    return perc, time_per_bin


def tuniv_ssfr_prior(objs, flds, pfile=None):
    # NEW TUNIV PRIOR
    ssfrs = []
    for i in range(len(objs)):
        model = pfile.load_model(objname=objs[i], field=flds[i])  # load_model prints zred, tuniv
        tuniv = WMAP9.age(model.params['zred']).value
        ssfrs.append(1 / tuniv)
        print(i)
    perc = np.percentile(ssfrs, [16., 50., 84.])
    print(perc)

    return perc


def draw_ssfr_from_prior(objs, flds, fig=None, axes=None, ndraw=1e4, alpha_sfh=1.0, pfile=None, show=True, t=None):
    # pfile=e_params or n_params

    # let's do it
    ndraw = int(ndraw)
    # zred = np.array([0.0, 0.5, 1.5, 2.5])  # where do we measure?
    zred = np.array([3.5])
    logmass = np.array([10.])
    smass_factor = 0.8  # add in a not-unreasonable stellar mass <--> total mass conversion
    # logmass *= 1.01  # note (10^(x*y) = 10^y / 0.8, for y in [9, 11] --> x ~= 1.01)
    minssfr, maxssfr = 1e-14, 1e-5

    fs = 20

    # for i, z in enumerate(zred):
    # for index, redshift in zred: i.e. for (0, 0.0), (1, 0.5), (2, 1.5), (3, 2.5)

    # new redshift, new model
    # model = pfile.load_model(zred=z, alpha_sfh=alpha_sfh, **pfile.run_params)  # load_model(objname, field):
    model = pfile.load_model(objname=objs[0], field=flds[0])
    agebins = model.params['agebins']

    # create prior & draw from prior
    prior = model._config_dict['sfr_fraction']['prior']
    print(prior)

    mass = np.zeros(shape=(agebins.shape[0], ndraw))
    '''
    for n in range(ndraw):
        mass[:, n] = pfile.sfrac_to_masses(logmass=logmass, sfrac=prior.sample(), agebins=agebins)
    '''
    # Jonathan-- this is the part you'll need to edit. i create a prior and sample from it with the lines above.
    # however you have a different setup than me now!
    # what you want from this calculation is an array of masses in each bin drawn directly from the prior
    # the two lines below will create a series of sfr_fractions drawn directly from the prior
    # what you have to do is convert these sfr_fractions in 'flatchain' into masses in each bin, in an array with
    # the same shape as the 'mass' vector above
    # these calculations are done in lines 525-533 here:
    # https://github.com/jhcohn/prospector-eelgs/blob/master/eelg_fixedmet_params.py
    # where each line in 'flatchain' below is a set of 'sfr_fractions'
    # Here you can assume the total log(M) = 10, as stated above
    flatchain = np.random.dirichlet(tuple(1.0 for x in xrange(6)), ndraw)
    # flatchain = (array([0.stuff, ..., 0.stuff]), ndraw)  # array is ndraw arrays of 6 different 0.stuffs
    flatchain = flatchain[:, :-1]

    # 525 - 533
    # use fractions()  # need to convert logmass to linear mass before doing line 533
    for n in range(ndraw):
        mass[:, n] = pfile.sfrac_to_masses(logmass=logmass, agebins=agebins, sfrac=flatchain[n])
        '''
        print(n)
        fractions = flatchain[n]  # np.array(model.params['sfr_fraction'])
        bin_fractions = np.append(fractions, (1 - np.sum(fractions)))
        time_per_bin = []
        for (t1, t2) in agebins:
            time_per_bin.append(10 ** t2 - 10 ** t1)
        # print(bin_fractions)
        # print(time_per_bin)
        bin_fractions *= np.array(time_per_bin)
        bin_fractions /= bin_fractions.sum()
        mass[:, n] = bin_fractions * pfile.transform_logmass_to_mass(logmass=model.params['mass'])
        '''
    # convert to sSFR
    # print(mass.shape)  # (6, ndraw)
    # print(mass)
    # sprint(mass[0, :].sum(axis=0))  # consistently ~5.4e12
    # print(agebins)  # [[0, 8], [8, 8.7],...]
    time_per_bin = np.diff(10**agebins, axis=-1)[:, 0]
    # print(np.diff(10**agebins, axis=-1))  # [[10^8 - 10^0],[10^8.7 - 10^8],...]
    # print(time_per_bin)  # [1e8, 4.01e8, 4.99e8, 2.84e8, 3.65e8, 4.685e8]
    # print(time_per_bin[0].sum)  # "<built-in method sum of numpy.float64 object at 0x7f99da4993d8>"
    # print(mass[0, :].sum(axis=0) / time_per_bin[0].sum() / 10**logmass)  # 5.3e-6 to 5.4e-6 > 1e-7 --> clipped
    # print(np.clip(mass[0, :].sum(axis=0) / time_per_bin[0].sum() / 10**logmass, minssfr, maxssfr))  ~ 5e-6 or 1e-7

    ssfr = np.zeros(shape=(len(mass), len(mass[0])))  # 6, ndraw
    # print(len(ssfr[0]))
    # print(mass[0, :].sum(axis=0))
    for i in range(len(mass)):  # 6
        # ssfr[i, :] = np.log10(mass[i, :] / time_per_bin[i] / 10**logmass)
        # ssfr[i, :] = np.log10(np.clip(mass[i, :].sum(axis=0) / time_per_bin[i].sum() / 10**logmass, minssfr, maxssfr))
        # ssfr[i, :] = np.log10(mass[i, :].sum(axis=0) / time_per_bin[i].sum() / (10**logmass / smass_factor))  # smass
        # ssfr[i, :] = np.log10(mass[i, :].sum(axis=0) / t / (10**logmass / smass_factor))  # smass
        ssfr[i, :] = np.log10(mass[i, :].sum(axis=0) / t / 10**logmass)  # smass

    # print(ssfr)
    print('perc', np.percentile(ssfr[0], [16, 50, 84]))
    # median of each bin should have ssfr of 1/tuniv at that redshift (i.e. continuous SFH)
    # print(ssfr)  always [-7.]
    # print(len(mass[0]), len(mass))  # 10000, 6

    '''
    ssfr_per_bin = []
    for i in range(6):
        ssfr_per_bin.append(mass[i, :].sum(axis=0) / time_per_bin[i].sum() / 10**logmass)
    # print(ssfr_per_bin)  # six ssfrs
    # print(np.average(agebins, axis=1))  # six times centered at each bin

    x = 10**np.average(agebins, axis=1) / 10**9
    x[0] = 3.5*10**-2
    axes.plot(x, ssfr_per_bin, 'o')
    axes.set_xlim(10**-2, 10**max(agebins[-1]) / 10**9)
    axes.set_ylim(10**-8, 10**-4)
    axes.set_yscale('log')
    axes.set_xscale('log')
    axes.set_ylabel(r'sSFR [yr$^{-1}$]', fontsize=20)
    axes.set_xlabel(r'Lookback time [Gyr]', fontsize=20)
    for listy in agebins:
        axes.axvline(x=10**listy[1] / 10**9, color='k', ls='-')
    plt.show()
    '''
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


def bootstrap(X, n=None, X_err=None):
    """
    Bootstrap resample an array_like
    Parameters
    :param X: array-like data to resample
    :param n: int, optional length of resampled array, equal to len(X) if n == None
    :param X_err:
    :return: X_resamples
    """
    if n is None:
        n = len(X)

    resample_i = np.floor(np.random.rand(n) * len(X)).astype(int)
    # creates len(n) array filled with random numbers from floor([0 to 1) * len(X))
    # i.e. it's a len(X) array filled with random integers from 0 to len(X)

    # resample_i=np.random.randint(low=0, high=len(X)-1, size=len(X))
    # print(resample_i)
    X_resample = X[resample_i]  # take X and use indices chosen randomly above
    if X_err != None:
        X_err_resample = X_err[resample_i]
        return X_resample, X_err_resample
    else:
        return X_resample, resample_i


def randraw(infile, logmass, num=1000):  # num=1000
    """
    For a given galaxy, randraw samples the posterior num times for each point in extra_output['extras']['sfh'][i]

    :param infile: ID_field_base_extra_out.py, where extra_output is stored using output.py
    :param num: number of times to sample the galaxy posterior at each point in the SFH
    :return: draw_from_sfh = 22 x num, lists the num random posterior samples, t_sfh = time vector associated with SFH
    """
    with open(infile, 'rb') as exout:
        extra_output = pickle.load(exout)

#    draw_from_sfh = np.zeros(shape=(len(extra_output['extras']['ssfr']), num))  # shape=(22, num)  # BUCKET MASS WRONG?
    draw_from_sfh = np.zeros(shape=(len(extra_output['extras']['sfh']), num))  # shape=(22, num)
    # print(len(extra_output['extras']['ssfr']), len(extra_output['extras']['ssfr'][0]))  # 22, 2000
    # print(len(draw_from_sfh), len(draw_from_sfh[0]))  # 22, num

    '''  # BUCKET MASS WRONG?
    for i in range(len(extra_output['extras']['ssfr'])):  # at each of these 22 points
        for j in range(num):  # randomly draw from the ssfr posterior num times
            draw_from_sfh[i][j] = extra_output['extras']['ssfr'][i][random.randint(0, num)] / 0.8  # BUCKET 1/0.8 here?
    '''
    for i in range(len(extra_output['extras']['sfh'])):  # at each of these 22 points
        for j in range(num):  # randomly draw from the ssfr posterior num times
            draw_from_sfh[i][j] = extra_output['extras']['sfh'][i][random.randint(0, num)] / (10**logmass) / 0.8

    return draw_from_sfh, extra_output['extras']['t_sfh']


'''
# Fine check, but not actual goal
def bootdraw(infile):
    """
    Attempting to combine randraw and bootstrap intelligently in order to use Leo's bootstrapping code

    :param infile: ID_field_base_extra_out.py, where extra_output is stored using output.py
    :return: draw_from_sfh = 22 x num, lists the num random posterior samples, t_sfh = time vector associated with SFH
    """
    with open(infile, 'rb') as exout:
        extra_output = pickle.load(exout)

    num = len(extra_output['extras']['ssfr'][0])
    # print(num, 'num')  # 2000
    draw_from_sfh = np.zeros(shape=(len(extra_output['extras']['ssfr']), num))  # shape=(22, num)

    for i in range(len(extra_output['extras']['ssfr'])):  # at each of these 22 points
        n = extra_output['extras']['ssfr'][i]
        # for j in range(len(n)):  # randomly draw from the ssfr posterior num times
        resample_i = np.floor(np.random.rand(len(n)) * len(n)).astype(int)
        draw_from_sfh[i] = n[resample_i]
        # draw_from_sfh[i][j] = extra_output['extras']['ssfr'][i][random.randint(0, num)]

    return draw_from_sfh, extra_output['extras']['t_sfh']
'''


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

    return all_draws


def plot_sfhs(draws, t, lw=1, elist=None, llist=None, uvj_in=False, spec=True, sigma=1, save=False, title=None,
              priors=None, show=False, tpbs=None, tuniv=False, nums=[1., 1.], both=False, seconds=None, ratio=False):
    """
    Plots SFH stacks for two different galaxy samples side-by-side

    :param percs: list of two smoothed percs, for two different galaxy samples, each output by smooth(perc)
    :param t: time vector output by randraw
    :param lw: line width
    :param spec: if stacking specific SFR instead of plain SFR, spec=True
    :param sigma: how many sigma of error we want to show in the plot
    :return: plot
    """
    log = 0
    x = [0.5e-8, 1e-8, 1.5e-8]  # , 2e-8]  # used if log=0
    # y = [0.0, 0.05, 0.10, 0.15]  # used regardless of log
    # y = [0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70]  # used regardless of log
    y = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.55]#, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]
    # , 0.35]  # used regardless of log; max 0.35 for 1e-9 bin spacing

    if ratio:
        fig = plt.figure()
        ax1 = plt.subplot(1, 1, 1)
    else:
        fig = plt.figure()
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2, sharey=ax1, sharex=ax1)
        ax1.set_xlim(xmin=3*10**-14, xmax=1.75e-8)  # 2.25e-8)  # 3*10**-12, 7e-8
        ax1.xaxis.set_ticks(x)
        ax1.set_ylim(ymin=0., ymax=0.57) #ymax=0.42)  # 0.34)  # 0.33)  # 0.28)  # 0.39)
        ax1.yaxis.set_ticks(y)
        fig.subplots_adjust(wspace=0)

    if not log:
        uvj_loc = 1  # 'upper right'
        labels = [5, 10, 15]  # , 20]
        ax1.set_xticklabels(labels)
    elif log:
        uvj_loc = 'upper center'
        ax1.set_xlim(xmin=3 * 10 ** -12, xmax=10 ** -7)  # 3*10**-12
        ax1.set_xscale("log")
        x = [1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]
        ax1.set_xlim(xmin=3*10**-15, xmax=1e-7)  # 3*10**-12
        ax1.xaxis.set_ticks(x)
    ax1.set_ylabel(r'Cumulative distribution', fontsize=30)

    if uvj_in:  # also requires elist, llist to be not None; this insets uvj plots onto the top right of plot!
        ht = 10*0.28  # 8*0.28
        wd = 10*0.32  # 8*0.32
        ecols = ['purple' for idx in range(len(elist))]
        lcols = ['blue' for idx in range(len(llist))]
        # for colr in range(len(elist)):
        #     cols.append('purple')
        # create inset axis: width (%), height (inches), location
        # loc=1 (upper right), loc=2 (upper left) --> loc=3 (lower left), loc=4 (lower right); loc=7 (center right)
        # https://stackoverflow.com/questions/10824156/matplotlib-legend-location-numbers
        if ratio:
            inset_axes(ax1, width=wd, height=ht, loc=10)  # 20%
            uvj.uvj_plot(-1, 'all', objlist=elist, title=False, labels=False, lims=True, size=20, show=False, col=ecols)
            inset_axes(ax1, width=wd, height=ht, loc=7)  # 20%
            uvj.uvj_plot(-1, 'all', objlist=llist, title=False, labels=False, lims=True, size=20, show=False,
                         col=lcols,
                         hist=True)
        else:
            inset_axes(ax1, width=wd, height=ht, loc=uvj_loc)  # 20%
            uvj.uvj_plot(-1, 'all', objlist=elist, title=False, labels=False, lims=True, size=20, show=False, col=ecols)
            inset_axes(ax2, width=wd, height=ht, loc=uvj_loc)  # 20%
            uvj.uvj_plot(-1, 'all', objlist=llist, title=False, labels=False, lims=True, size=20, show=False, col=lcols,
                         hist=True)

    if priors is not None:
        if not ratio:
            ax1.axvline(x=priors[0][1], color='k', ls='--', label='Prior')  # median, +/- 1sigma for EELG prior
            ax2.axvline(x=priors[1][1], color='k', ls='--', label='Prior')  # median, +/- 1sigma for SFG prior

    # DO THE PLOTTING! AND KS TEST!
    print('hi')
    logged1 = []
    logged2 = []
    for i in range(len(draws[0])):
        logged1.append(np.log10(draws[0][i]))
    for i in range(len(draws[1])):
        logged2.append(np.log10(draws[1][i]))
    perc1 = np.percentile(draws[0], [25, 50, 75])
    perc2 = np.percentile(draws[1], [25, 50, 75])
    # print(perc1, perc2)

    # KS test bins: 100 bins corresponds to bin width 1e-9 (1 Gyr^-1)
    num_bins = 100
    hi = np.histogram(draws[0], bins=num_bins, range=(10**-12, 10**-7))
    hey = np.histogram(draws[1], bins=num_bins, range=(10**-12, 10**-7))
    print(hi[0], hey[0])
    h1 = [x1*127*4/3/19 for x1 in hi[0]]
    e1 = [x2*4/3 for x2 in hey[0]]
    print('KS', scipy.stats.ks_2samp(h1, e1))  # ~0.005
    print('A-D', scipy.stats.anderson_ksamp((h1, e1)))  # ~5e-5

    hi2 = np.histogram(seconds[0], bins=num_bins, range=(10 ** -12, 10 ** -7))
    hey2 = np.histogram(seconds[1], bins=num_bins, range=(10 ** -12, 10 ** -7))
    print(hi2[0], hey2[0])
    h2 = [y1*127/19 for y1 in hi[0]]
    e2 = [y2 for y2 in hey[0]]
    print('KS2', scipy.stats.ks_2samp(h2, e2))  # 0.01 (really?? maybe bc 4e3 pts, instead of just 3e3 above?)
    print('A-D2', scipy.stats.anderson_ksamp((hi2[0], hey2[0])))  # 1.76e-5

    # NOTE: the above numbers changed WRONG: NOW CORRECT, for anderson-darling: (100: 0.043)
    # print(min(draws[1]), max(draws[1]))  # ~4.7e-14, 21*1e-9
    print(np.percentile(draws[1], [16, 50, 84]))  # 4.08e-10, 1.38e-9, 2.89e-9
    # print(min(draws[0]), max(draws[0]))  # 2.5e-10, 12e-9
    print(np.percentile(draws[0], [16, 50, 84]))  # 1.25e-9, 3.05e-9, 4.59e-9
    # sqrt(N_SFG) = sqrt(167) = 12.92285 --> error on mean: {1.38e-9}^{+0.117e-9}_{-0.00752e-9} (max 1.497)
    # sqrt(N_EELG) = sqrt(18) = 4.24264 --> error on mean: {3.05e-9}^{+0.363e-9}_{-0.424e-9} (min 2.626)

    display_num = np.arange(1e-14, 1e-7 + 1e-9, 0.67e-9)  # display_num = np.arange(1e-11, 1e-7 + 1e-9, 1e-9)
    # use 1e-9 OR 0.67e-9 for bin width in np.arange^
    # display_num = np.linspace(1e-12, 1e-7, num=30)
    fs = 20
    if not log:
        if both:
            # print(len(draws[1]), len(seconds[1]))  # = len(SFGs) * 1000
            hi3 = np.histogram(draws[1], bins=num_bins, range=(10 ** -12, 10 ** -7))
            hey3 = np.histogram(seconds[1], bins=num_bins, range=(10 ** -12, 10 ** -7))
            print('KS3', scipy.stats.ks_2samp(hi3[0], hey3[0]))  # ~0.19 (believable)
            print('A-D3', scipy.stats.anderson_ksamp((hi3[0], hey3[0])))  # 0.17 <-- NOT SIGNIFICANT

            # print(len(draws[0]), len(seconds[0]))  # = len(EELGs) * 1000
            hi4 = np.histogram(draws[0], bins=num_bins, range=(10 ** -12, 10 ** -7))
            hey4 = np.histogram(seconds[0], bins=num_bins, range=(10 ** -12, 10 ** -7))
            print('KS4', scipy.stats.ks_2samp(hi4[0], hey4[0]))  # ~0.9995 (SHOULD BE VERY DIFF)
            print('A-D4', scipy.stats.anderson_ksamp((hi4[0], hey4[0])))  # ~0.76 WTF THIS IS NOT CORRECT
            # ax1.set_ylabel(r'Cumulative distribution', fontsize=30)
            ax1.set_ylabel(r'Fraction of composite SED', fontsize=30)
            # y = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # , 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]
            # ax1.set_ylim(ymin=0., ymax=1.0)  # ymax=0.42)  # 0.34)  # 0.33)  # 0.28)  # 0.39)
            # ax1.set_xlim(xmin=3 * 10 ** -14, xmax=6e-8)  # 2.25e-8)  # 3*10**-12, 7e-8
            # ax1.xaxis.set_ticks(x)

            ax1.hist(draws[0], histtype='bar', bins=display_num, weights=[1. / (len(draws[0]))] * len(draws[0]),
                     color='purple', alpha=0.5, lw=2, label='EELGs, most recent bin')
            ax2.hist(draws[1], histtype='bar', bins=display_num, weights=[1. / (len(draws[1]))] * len(draws[1]),
                     color='b', alpha=0.5, lw=2, label='SFGs, most recent bin')
            ax1.hist(seconds[0], histtype='bar', bins=display_num, weights=[1. / (len(seconds[0]))] * len(seconds[0]),
                     color='purple', ec='purple', hatch='/', fill=False, lw=2, label='EELGs, second bin')  # ls='--',
            ax2.hist(seconds[1], histtype='bar', bins=display_num, weights=[1. / (len(seconds[1]))] * len(seconds[1]),
                     color='b', ec='b', hatch='/', lw=2, fill=False, label='SFGs, second bin')  # ls='--',

            fig.text(0.5, 0.04, r'sSFR [Gyr$^{-1}$]', ha='center', fontsize=30)  # 30
            ax1.legend(numpoints=1, loc='lower left', bbox_to_anchor=(0.05, 0.842), prop={'size': fs})  # (0.05, 0.88)
            ax2.legend(numpoints=1, loc='lower left', bbox_to_anchor=(0.05, 0.842), prop={'size': fs})  # (0.05, 0.88)
            # ax1.legend(numpoints=1, loc='upper left', prop={'size': fs})  # (0.05, 0.88)
            # ax2.legend(numpoints=1, loc='upper left', prop={'size': fs})  # (0.05, 0.88)  # bucket just now

            fs_ticks = 25
            # ax1.set_xticklabels([r'$10$', r'$20$', r'$30$',  r'$40$', r'$50$'], size=fs_ticks)
            ax2.set_xticklabels([r'$5$', r'$10$', r'$15$'], size=fs_ticks)  # bucket just now
            ax2.set_xticklabels([r'$5$', r'$10$', r'$15$'], size=fs_ticks)  # bucket just now
            ax1.set_yticklabels([r'$0.00$', r'$0.05$', r'$0.10$', r'$0.15$', r'$0.20$', r'$0.25$', r'$0.30$',
                                 r'$0.35$', r'$0.40$', r'$0.45$', r'$0.50$', r'$0.55$'], size=fs_ticks)
            # ax1.set_yticklabels([r'$0.00$', r'$0.10$', r'$0.20$', r'$0.30$', r'$0.40$', r'$0.50$', r'$0.60$',
            #                      r'$0.70$', r'$0.80$', r'$0.90$', r'$1.0$'], size=fs_ticks)

        elif ratio:
            x = [0., 2., 4., 6., 8., 10., 12., 14., 16., 18.]  # , 2e-8]  # used if log=0
            #y = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]  # bucket just now
            y = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # , 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]
            ax1.set_xlim(xmin=0., xmax=19.)  # 2.25e-8)  # 3*10**-12, 7e-8
            ax1.xaxis.set_ticks(x)
            ax1.set_ylim(ymin=0., ymax=1.0) #0.37  # ymax=0.42)  # 0.34)  # 0.33)  # 0.28)  # 0.39)
            ax1.yaxis.set_ticks(y)
            # print(np.percentile(np.histogram(seconds[0], bins=np.arange(0, 20., 0.5)), [16., 50., 84.]))
            print(max(seconds[0]), max(seconds[1]), 'maxes')
            ax1.hist(seconds[0], histtype='step', bins=np.arange(0, 20., 0.5), cumulative=True, normed=1,
                     color='purple', alpha=0.75, lw=2, label='EELGs')
            # weights=[1. / (nums[0] * 10 ** 3)] * len(seconds[0]),
            ax1.hist(seconds[1], histtype='step', bins=np.arange(0, 20., 0.5), cumulative=True, normed=1,  # bucket just now
                     color='b', alpha=0.75, lw=2, label='SFGs')
            # weights=[1. / (nums[1] * 10 ** 3)] * len(seconds[1]),
            fig.text(0.5, 0.04, r'(sSFR most recent bin) / (sSFR, second bin)', ha='center', fontsize=30)  # 30
            fs_ticks = 25
            ax1.set_xticklabels([r'$0$', r'$2$', r'$4$', r'$6$', r'$8$', r'$10$', r'$12$', r'$14$', r'$16$', r'$18$'],
                                size=fs_ticks)
            #ax2.set_xticklabels([r'$0$', r'$2$', r'$4$', r'$6$', r'$8$', r'$10$', r'$12$', r'$14$', r'$16$', r'$18$'],
            #                    size=fs_ticks)  # bucket just now
            #ax1.set_yticklabels([r'$0.00$', r'$0.05$', r'$0.10$', r'$0.15$', r'$0.20$', r'$0.25$', r'$0.30$',
            #                     r'$0.35$'], size=fs_ticks)  # bucket just now
            ax1.set_yticklabels([r'$0.00$', r'$0.10$', r'$0.20$', r'$0.30$', r'$0.40$', r'$0.50$', r'$0.60$',
                                 r'$0.70$', r'$0.80$', r'$0.90$', r'$1.0$'], size=fs_ticks)
            ax1.axvline(x=1., color='k', ls='--')
            #ax2.axvline(x=1., color='k', ls='--')  # bucket just now
            ax1.legend(numpoints=1, loc='lower left', bbox_to_anchor=(0.8, 0.742), prop={'size': fs})  # 0.88
            #ax2.legend(numpoints=1, loc='lower left', bbox_to_anchor=(0.75, 0.842), prop={'size': fs})  # 0.88  # bucket just now
        else:
            ax1.hist(draws[0], histtype='bar', bins=display_num, weights=[1. / (nums[0]*len(draws[0]))] * len(draws[0]),
                     color='purple', alpha=0.75, lw=2, label='EELGs')
            ax2.hist(draws[1], histtype='bar', bins=display_num, weights=[1. / (nums[1]*len(draws[1]))] * len(draws[1]),
                     color='b', alpha=0.75, lw=2, label='SFGs')
            fig.text(0.5, 0.04, r'sSFR (most recent bin) [Gyr$^{-1}$]', ha='center', fontsize=30)  # 30
            ax1.legend(numpoints=1, loc='lower left', bbox_to_anchor=(0.05, 0.88), prop={'size': fs})
            ax2.legend(numpoints=1, loc='lower left', bbox_to_anchor=(0.05, 0.88), prop={'size': fs})
            fs_ticks = 25
            ax1.set_xticklabels([r'$5$', r'$10$', r'$15$'], size=fs_ticks)
            ax2.set_xticklabels([r'$5$', r'$10$', r'$15$'], size=fs_ticks)
            ax1.set_yticklabels([r'$0.00$', r'$0.05$', r'$0.10$', r'$0.15$', r'$0.20$', r'$0.25$', r'$0.30$', r'$0.35$',
                                 r'$0.40$', r'$0.45$', r'$0.50$', r'$0.55$'], size=fs_ticks)

        # ax1.legend(numpoints=1, loc='lower left', bbox_to_anchor=(0.05, 0.842), prop={'size': fs})  # (0.05, 0.88)#
        # ax2.legend(numpoints=1, loc='lower left', bbox_to_anchor=(0.05, 0.842), prop={'size': fs})  # (0.05, 0.88)
        # ax1.legend(numpoints=1, loc='upper left', prop={'size': fs})
        # ax2.legend(numpoints=1, loc='upper left', prop={'size': fs})
    else:
        ax1.hist(draws[0], histtype='bar', bins=display_num, weights=[3./(18*3*10**3)]*len(draws[0]), color='b',
                 alpha=0.5, lw=2, label='EELGs')
        ax2.hist(draws[1], histtype='bar', bins=display_num, weights=[1./(87*3*10**3)]*len(draws[1]), color='r',
                 alpha=0.5, lw=2, label='SFGs')
        fig.text(0.5, 0.04, 'SSFR (most recent bin; yr$^{-1}$)', ha='center', fontsize=30)  # 30
        ax1.legend(numpoints=1, loc='upper right', prop={'size': fs})
        ax2.legend(numpoints=1, loc='upper right', prop={'size': fs})

    plt.setp(ax2.get_yticklabels(), visible=False)  # hide y-axis labels on right-hand subplot to prevent overlap # bucket just now
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size': 22})

    plt.show()


# Define a function to make the ellipses
def ellipse(ra, rb, ang, x0, y0, Nb=100):
    xpos, ypos = x0, y0
    radm, radn = ra, rb
    an = ang
    co, si = np.cos(an), np.sin(an)
    the = np.linspace(0, 2 * np.pi, Nb)
    X = radm * np.cos(the) * co - si * radn * np.sin(the) + xpos
    Y = radm * np.cos(the) * si + co * radn * np.sin(the) + ypos
    return X, Y


def percellipse(percs_a, percs_b, ang, Nb=100):
    xpos, ypos = percs_a[1], percs_b[1]
    radm, radn = percs_a[2] - percs_a[0], percs_b[2] - percs_b[0]
    an = ang
    co, si = np.cos(an), np.sin(an)
    the = np.linspace(0, 2 * np.pi, Nb)
    X = radm * np.cos(the) * co - si * radn * np.sin(the) + xpos
    Y = radm * np.cos(the) * si + co * radn * np.sin(the) + ypos
    return X, Y


def twodplot(recent_x, second_y, histcol):
    #plt.ion()

    # Define the x and y data
    # For example just using random numbers
    x = recent_x
    y = second_y

    # Set up default x and y limits
    xlims = [min(x), max(x)]
    ylims = [min(y), max(y)]
    # xlims = [min([0., min(x)]), min([20., max(x)])]
    # ylims = [min([0., min(y)]), min([20., max(y)])]
    # xlims = [0., 20.]
    # ylims = [0., 20.]

    # Set up your x and y labels
    xlabel = r'SSFR, most recent bin [Gyr$^{-1}$]'
    ylabel = r'SSFR, second most recent bin [Gyr$^{-1}$]'

    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left + width + 0.02

    # Set up the geometry of the three plots
    rect_temperature = [left, bottom, width, height]  # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25]  # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height]  # dimensions of y-histogram

    # Set up the size of the figure
    fig = plt.figure(1, figsize=(9.5, 9))

    # Make the three plots
    axTemperature = plt.axes(rect_temperature)  # temperature plot
    axHistx = plt.axes(rect_histx)  # x histogram
    axHisty = plt.axes(rect_histy)  # y histogram

    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Find the min/max of the data
    xmin = min(xlims)
    xmax = max(xlims)
    ymin = min(ylims)
    ymax = max(y)

    # Make the 'main' temperature plot
    # Define the number of bins
    nxbins = 50
    nybins = 50
    nbins = 100

    xbins = np.linspace(start=xmin, stop=xmax, num=nxbins)
    ybins = np.linspace(start=ymin, stop=ymax, num=nybins)
    xcenter = (xbins[0:-1] + xbins[1:]) / 2.0
    ycenter = (ybins[0:-1] + ybins[1:]) / 2.0
    aspectratio = 1.0 * (xmax - 0) / (1.0 * ymax - 0)

    H, xedges, yedges = np.histogram2d(y, x, bins=(ybins, xbins))
    X = xcenter
    Y = ycenter
    Z = H

    # Plot the temperature data
    cax = (axTemperature.imshow(H, extent=[xmin, xmax, ymin, ymax],
                                interpolation='nearest', origin='lower', aspect=aspectratio))

    # Plot the temperature plot contours
    contourcolor = 'white'
    # '''
    xcenter = np.mean(x)
    ycenter = np.mean(y)
    ra = np.std(x)
    rb = np.std(y)
    '''
    xcenter = np.median(x)
    ycenter = np.median(y)
    ra = (np.percentile(x, [84.]) - np.percentile(x, [16.])) / 2
    rb = (np.percentile(y, [84.]) - np.percentile(y, [16.])) / 2
    '''
    ang = 0

    X, Y = ellipse(ra, rb, ang, xcenter, ycenter)
    axTemperature.plot(X, Y, "w:", ms=1, linewidth=2.0)
    axTemperature.annotate('$1\\sigma$', xy=(X[15], Y[15]), xycoords='data', xytext=(10, 10),
                           textcoords='offset points', horizontalalignment='right',
                           verticalalignment='bottom', fontsize=25, color=contourcolor)

    X, Y = ellipse(2 * ra, 2 * rb, ang, xcenter, ycenter)
    axTemperature.plot(X, Y, "w:", color=contourcolor, ms=1, linewidth=2.0)
    axTemperature.annotate('$2\\sigma$', xy=(X[15], Y[15]), xycoords='data', xytext=(10, 10),
                           textcoords='offset points', horizontalalignment='right',
                           verticalalignment='bottom', fontsize=25, color=contourcolor)

    X, Y = ellipse(3 * ra, 3 * rb, ang, xcenter, ycenter)
    axTemperature.plot(X, Y, "w:", color=contourcolor, ms=1, linewidth=2.0)
    axTemperature.annotate('$3\\sigma$', xy=(X[15], Y[15]), xycoords='data', xytext=(10, 10),
                           textcoords='offset points', horizontalalignment='right',
                           verticalalignment='bottom', fontsize=25, color=contourcolor)

    # Plot the axes labels
    axTemperature.set_xlabel(xlabel, fontsize=25)
    axTemperature.set_ylabel(ylabel, fontsize=25)

    # Make the tickmarks pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family('serif')

    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family('serif')

    # Set up the plot limits
    axTemperature.set_xlim(xlims)
    axTemperature.set_ylim(ylims)

    # Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax - xmin) / nbins)
    ybins = np.arange(ymin, ymax, (ymax - ymin) / nbins)

    # Plot the histograms
    axHistx.hist(x, bins=xbins, color=histcol)
    axHisty.hist(y, bins=ybins, orientation='horizontal', color=histcol)

    # Set up the histogram limits
    axHistx.set_xlim(min(x), max(x))
    axHisty.set_ylim(min(y), max(y))

    # Make the tickmarks pretty
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')

    # Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')

    # Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))

    # Show the plot
    plt.show()


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def simpler(recentx, secondy, col):
    # the random data
    # x = np.random.randn(1000)
    # y = np.random.randn(1000)
    x = recentx[0]
    y = secondy[0]
    x2 = recentx[1]
    y2 = secondy[1]

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    # added
    left2, width2 = 0.55, 0.4
    bottom2, height2 = 0.55, 0.4
    # bottom_h2 = left_h2 = left2 + width2 + 0.02
    rect_scatter2 = [left2, bottom2, width2, height2]

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    fig1 = plt.figure(1, figsize=(8, 8))  # right now I have to manually set aspect by manually expanding the image...

    #axScatter = plt.axes(rect_scatter)
    #axHistx = plt.axes(rect_histx)
    #axHisty = plt.axes(rect_histy)
    specs = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[3, 1], height_ratios=[1, 3], wspace=0.05, hspace=0.05)

    axScatter = fig1.add_subplot(specs[1, 0])
    axHistx = fig1.add_subplot(specs[0, 0])
    axHisty = fig1.add_subplot(specs[1, 1])

    # axScatter.set_aspect('equal')
    # axHistx.set_aspect('equal')
    # axHisty.set_aspect('equal')

    # added
    # axScatter2 = plt.axes(rect_scatter2)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:

    # now determine nice limits by hand:
    binwidth = 0.25
    xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
    lim = (int(xymax / binwidth) + 1) * binwidth

    xlims = [0., 17.]
    ylims = [0., 17.]

    # axScatter.set_xlim((-lim, lim))
    # axScatter.set_ylim((-lim, lim))
    # axScatter.set_xlim(0., max([max(x), max(x2)]))
    # axScatter.set_ylim(0., max([max(y), max(y2)]))
    axScatter.set_xlim(xlims[0], xlims[1])
    axScatter.set_ylim(ylims[0], ylims[1])

    # axScatter2.set_xlim(0., max([max(x), max(x2)]))
    # axScatter2.set_ylim(0., max([max(y), max(y2)]))

    '''
    axScatter.scatter(x2, y2, alpha=0.01, color=col[1])
    axScatter.scatter(x, y, alpha=0.01, color=col[0])
    '''
    # OR INSTEAD OF SCATTER:
    '''
    import ewmass as ewm
    # H, xedges, yedges = np.histogram2d(y, x, bins=(ybins, xbins))
    # H2, xedges2, yedges2 = np.histogram2d(y, x, bins=(ybins, xbins))
    ewm.plotter(x2, y2, color='b', xs=[0., 20.], ys=[0., 20.], ax=axScatter, scat=False,
                levels=[1.-0.9973, 1-0.9545, 1-0.6827, 1.], alpha=0.9)
    ewm.plotter(x, y, color='purple', xs=[0., 20.], ys=[0., 20.], ax=axScatter, scat=False,
                levels=[1.-0.9973, 1-0.9545, 1-0.6827, 1.], alpha=0.9)
    '''
    cmap2 = plt.get_cmap('Blues')
    new_cmap2 = truncate_colormap(cmap2, 0.3, 2.0)
    axScatter.hist2d(x2, y2, bins=100, cmap=new_cmap2, cmin=len(x2)/1000)  # bins with counts < cmin won't be displayed!
    cmap = plt.get_cmap('Purples')
    new_cmap = truncate_colormap(cmap, 0.3, 2.0)
    axScatter.hist2d(x, y, bins=100, cmap=new_cmap, cmin=len(x)/1000)  # bins with counts < cmin won't be displayed!
    axScatter.plot(xlims, ylims, ls='--', color='k')
    # axScatter.gca().set_aspect('equal', adjustable='box')
    # print(len(x), len(x2), 'lens!')  # 19000, 127000

    bins = np.arange(-lim, lim + binwidth, binwidth)
    axHistx.hist(x, bins=bins, color=col[0], alpha=0.3, weights=[1./len(x)] * len(x))
    axHisty.hist(y, bins=bins, orientation='horizontal', color=col[0], alpha=0.3, weights=[1./len(y)] * len(y))
    axHistx.hist(x2, bins=bins, hatch='/', ec=col[1], fill=False, weights=[1./len(x2)] * len(x2))
    axHisty.hist(y2, bins=bins, orientation='horizontal', hatch='/', ec=col[1], fill=False, weights=[1./len(y2)] *
                 len(y2))
    axScatter.set_xlim(xlims[0], xlims[1])
    axScatter.set_ylim(ylims[0], ylims[1])
    axHistx.set_xlim(axScatter.get_xlim())
    axHistx.set_ylim(0., 0.25)
    axHisty.set_ylim(axScatter.get_ylim())
    axHisty.set_xlim(0., 0.25)

    # Set up your x and y labels
    xlabel = r'<SSFR>$_{0-50}$n [Gyr$^{-1}$]'  # r'SSFR, most recent bin [Gyr$^{-1}$]'
    ylabel = r'<SSFR>$_{50-100}$n [Gyr$^{-1}$]'  # r'SSFR, second most recent bin [Gyr$^{-1}$]'
    axScatter.set_xlabel(xlabel, fontsize=30)
    axScatter.set_ylabel(ylabel, fontsize=30)

    # Make the tickmarks pretty
    fs = 20
    fs_ticks = 25
    xtick = [0., 5., 10., 15.]  # , 2e-8]  # used if log=0
    ytick = [0., 5., 10., 15.]  # , 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]
    axScatter.xaxis.set_ticks(xtick)
    axScatter.yaxis.set_ticks(ytick)
    axScatter.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    axScatter.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    axScatter.set_xticklabels([r'$0$', r'$5$', r'$10$', r'$15$'], size=fs_ticks)
    axScatter.set_yticklabels([r'$0$', r'$5$', r'$10$', r'$15$'], size=fs_ticks)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)

    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fs_ticks)
        label.set_family('serif')

    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fs_ticks)
        label.set_family('serif')

    # Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))

    plt.show()


if __name__ == "__main__":

    boot = 0

    corr = 0
    fico = 1
    newsfg = 0

    if corr:
        base = ['corr', 'corr']
        folders = ['pkl_ecorr/', 'pkl_ncorr/']
        mass = [9.48, 10.12]
        import eelg_varymet_params as e_params
        import eelg_varymet_params as n_params
        normal = True
    elif fico:
        base = ['fico', 'fico']
        folders = ['pkl_efico/', 'pkl_nfico/']
        mass = [9.42, 10.13]
        import eelg_fifty_params as e_params
        import eelg_fifty_params as n_params
        normal = True
    elif newsfg:
        base = ['fico', 'newsfg']
        folders = ['pkl_efico/', 'pkl_nnewsfg/']
        mass = [9.42, 9.84]
        import eelg_fifty_params as e_params
        import eelg_fifty_params as n_params
        normal = False

    import stellar_ages as sa
    e_objs, e_fields, l_objs, l_fields = sa.get_gal_lists(base, objlists=True, normal=normal)
    eelgs, lbgs = sa.get_gal_lists(base, objlists=False, normal=normal)

    pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[0]
    l_pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[1]

    # START STACKING
    t1 = []
    draws = []
    boots = []
    nummy = 0
    c = 0

    draws2 = []
    numl = 0
    cl = 0
    # t2 = []
    for glxy in eelgs:
        c += 1
        # file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
        file = pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            nummy += 1
            # temp = randraw(file)  # temp[0] lists num=1000 random posterior samples; temp[1] = time vector
            temp = randraw(file, mass[0])  # temp[0] lists num=1000 random posterior samples; temp[1] = time vector
            # temp = bootdraw(file)  # temp[0] lists num=1000 random posterior samples; temp[1] = time vector
            draws.append(temp[0])  # append random draw of ssfr
            # boots.append(bootstrap(temp[0]))
            t1.append(temp[1])
            print(len(temp[0]))  # 22, 1000
        else:
            print(file)

    lfile_loc = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/sfg_ssfrs'
    for glxy in lbgs:
        cl += 1
        # file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
        file = l_pkls + glxy + '_extra_out.pkl'

        if os.path.exists(file):
            numl += 1
            # temp = randraw(file)
            temp = randraw(file, mass[1])
            # temp = bootdraw(file)
            draws2.append(temp[0])
            # t2.append(temp[1])
        else:
            print(file)

    # stacker2(draws, t1)
    sig = 1  # what sigma error to show on plot
    all1 = stacker(draws, sigma=sig)
    all2 = stacker(draws2, sigma=sig)

    means1 = np.zeros(shape=(10**3))
    means2 = np.zeros(shape=(10**3))
    # len(gal_draws) = number of galaxies in stack; len(gal_draws[0]) = 22, len(gal_draws[0][0]) = num (1000)
    # print(len(draws), len(draws[0]), len(draws[0][0]))  # 19, 22, 1000
    for nu in range(1000):
        each = []
        each2 = []
        for gal in range(len(draws)):
            each.append(draws[gal][0][nu])
        for gal2 in range(len(draws2)):
            each2.append(draws2[gal2][0][nu])
        means1[nu] = np.mean(each)
        means2[nu] = np.mean(each2)
    print('errors on mean:')
    print(np.percentile(means1, [16., 50., 84.]))
    print(np.percentile(means2, [16., 50., 84.]))

    rec1 = []
    sec1 = []
    for i in range(len(all1[0])):
        recent = []
        second = []
        for j in (0, 1, 2):
            recent.append(all1[j][i])
        for k in (3, 4, 5, 6):
            second.append(all1[k][i])
        rec1.append(sum(recent) * 10**9 / 3)
        sec1.append(sum(second) * 10**9 / 4)
    rec2 = []
    sec2 = []
    for i in range(len(all2[0])):
        recent = []
        second = []
        for j in (0, 1, 2):
            recent.append(all2[j][i])
        for k in (3, 4, 5, 6):
            second.append(all2[k][i])
        rec2.append(sum(recent) * 10**9 / 3)
        sec2.append(sum(second) * 10**9 / 4)

    randome1 = []
    randome2 = []
    idx = 4
    for i in range(len(draws[0][0])):
        recent = []
        for j in (0, 1, 2):
            recent.append(draws[idx][j][i])
        randome1.append(sum(recent) * 10**9 / 3)
    for i in range(len(draws[0][0])):
        second = []
        for k in (3, 4, 5, 6):
            second.append(draws[idx][k][i])
        randome2.append(sum(second) * 10**9 / 4)

    randoms1 = []
    randoms2 = []
    idx2 = 0
    for i in range(len(draws2[0][0])):
        recent = []
        for j in (0, 1, 2):
            recent.append(draws2[idx2][j][i])
        randoms1.append(sum(recent) * 10**9 / 3)
    for i in range(len(draws2[0][0])):
        second = []
        for k in (3, 4, 5, 6):
            second.append(draws2[idx2][k][i])
        randoms2.append(sum(second) * 10**9 / 4)

    # twodplot(randoms1, randoms2, 'b')
    # twodplot(randome1, randome2, 'purple')

    # smooth_percs = perc1, perc2
    print(nummy, c, 'nume')
    print(numl, cl, 'numl')
    # smooth_percs = [smooth(perc1), smooth(perc2)]

    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)

    recents = [rec1, rec2]
    seconds = [sec1, sec2]

    simpler(recents, seconds, col=['purple', 'b'])
    # simpler(rec2, sec2, col='b')

    twodplot(rec1, sec1, 'purple')
    twodplot(rec2, sec2, 'b')
