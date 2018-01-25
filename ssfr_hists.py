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

'''
# Also plots it
def stacker2(gal_draws, times, sigma=1, spec=True, lw=1):
    """
    stacker2 takes input of random points drawn from a list of galaxies' SFH posteriors, concatenates them within each
    bin, and then calculates the median and 1 sigma errors in each bin
    gal_draws should be in format draws = [draw_from_sfh1, draw_from_sfh2, ...]
    times should be in format times = [t1, t2, t3, ...]
    each draw_from_sfh has shape=(22,num)

    :param gal_draws: list comprised of draw_from_sfh (i.e. comprised of the output from randraw) for a list of galaxies
    :param times: list comprised of extra_output['extras']['t_sfh'] that correspond respectively to entries in gal_draws
    :param sigma: how many sigma of error we want to show in the plot
    :return: perc = stored lists of the median and +/- sigma SFH values calculated from the gal_draws
    """

    # print(len(perc), len(perc[0]), len(perc[0][0]))  # number of galaxies (12), number of points (22), 2*sigma+1 (3)
    # print(len(gal_draws[0][0]))  # 1000

    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)

    for k in range(len(gal_draws)):
        perc = np.zeros(shape=(len(gal_draws[0]), 2 * sigma + 1))  # len(gal_draws[0])=22=len(t); len(perc)=22, len(perc[0])=3
        # append the num=1000 values in each gal_draws[k] at each of the 22 points to all_draws:
        for jj in xrange(len(times[k])):
            print(perc[jj, :])
            perc[jj, :] = np.percentile(gal_draws[k][jj, :], [16.0, 50.0, 84.0])

        ax1.plot(times[k], perc[:, 1], '-', color='k', lw=lw)  # median
        ax1.fill_between(times[k], perc[:, 0], perc[:, 2], color='k', alpha=0.3)  # fill region between +/- 1sigma
        ax1.plot(times[k], perc[:, 0], '-', color='k', alpha=0.3, lw=lw)  # -1sigma
        ax1.plot(times[k], perc[:, 2], '-', color='k', alpha=0.3, lw=lw)  # +1sigma

    ax1.set_yscale("log")
    ax1.set_xscale("log")
    plt.show()

    return perc
'''


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
              priors=None, show=False, tpbs=None, tuniv=False):
    """
    Plots SFH stacks for two different galaxy samples side-by-side

    :param percs: list of two smoothed percs, for two different galaxy samples, each output by smooth(perc)
    :param t: time vector output by randraw
    :param lw: line width
    :param spec: if stacking specific SFR instead of plain SFR, spec=True
    :param sigma: how many sigma of error we want to show in the plot
    :return: plot
    """
    label = r'Fraction'  # r'Stacked SFR$_{bin}$ / M$_{tot}$ [yr$^{-1}$]'
    x = [10**-8, 5*10**-8, 10**-7]  # used if log=0
    log = 1

    fig = plt.figure()
    ax1 = plt.subplot(1, 2, 1)
    ax2 = plt.subplot(1, 2, 2, sharey=ax1, sharex=ax1)
    if log:
        ax1.set_xscale("log")
        x = [10**-11, 10**-10, 10**-9, 10**-8, 10**-7]
    ax1.set_xlim(xmin=3*10**-12, xmax=10**-7)
    ax1.set_ylim(ymin=0., ymax=0.40)
    ax1.xaxis.set_ticks(x)
    ax1.set_ylabel(label)
    fig.subplots_adjust(wspace=0)
    # ax2.set_xscale("log")
    # ax2.set_xlim(xmin=10**-12, xmax=10**-7)
    # ax2.set_ylim(ymin=0., ymax=1.)

    if uvj_in:  # also requires elist, llist to be not None; this insets uvj plots onto the top right of plot!
        inset_axes(ax1, width=8*0.32, height=8*0.28, loc=1)  # 20%
        uvj.uvj_plot(-1, 'all', objlist=elist, title=False, labels=False, lims=True, size=20, show=False)
        # create inset axis: width (%), height (inches), location
        # loc=1 (upper right), loc=2 (upper left) --> loc=3 (lower left), loc=4 (lower right); loc=7 (center right)
        # https://stackoverflow.com/questions/10824156/matplotlib-legend-location-numbers
        inset_axes(ax2, width=8*0.32, height=8*0.28, loc=1)  # 20%
        uvj.uvj_plot(-1, 'all', objlist=llist, title=False, labels=False, lims=True, size=20, show=False)

    # DO THE PLOTTING!
    print('hi')
    logged1 = []
    logged2 = []
    for i in range(len(draws[0])):
        logged1.append(np.log10(draws[0][i]))
    for i in range(len(draws[1])):
        logged2.append(np.log10(draws[1][i]))
    perc1 = np.percentile(draws[0], [25, 50, 75])
    perc2 = np.percentile(draws[1], [25, 50, 75])
    print(perc1, perc2)

    hi = np.histogram(draws[0], bins=39, range=(10**-12, 10**-7))  # get p <0.05 at ~35 bins
    hey = np.histogram(draws[1], bins=39, range=(10**-12, 10**-7))

    print(scipy.stats.ks_2samp(hi[0], hey[0]))
    ax1.hist(draws[0], histtype='bar', bins=10, weights=[3./(18*3*10**3)]*len(draws[0]), color='b', alpha=0.5, lw=2,
             label='EELGs')
    ax2.hist(draws[1], histtype='bar', bins=10, weights=[1./(87*3*10**3)]*len(draws[1]), color='r', alpha=0.5, lw=2,
             label='SFGs')

    # ax1.set_yscale("log")
    '''
    ax1.set_ylim(ymin, 10**-7)  # , ymax
    ax1.set_xlim(10**-2, 2.5)  # (0, 2.5)  # (10**-2, 13.6)
    ax1.set_ylabel(label, fontsize=30)  # 30
    ax1.text(2*10**-2, 4*10**-8, 'EELGs', fontsize=30)  # use if uvj_in
    # ax1.text(0.5, 5*10**-8, 'EELGs', fontsize=30)  # use if uvj_in
    '''
    # ax2.set_yscale("log")
    '''
    ax2.set_ylim(ymin, 10**-7)  # , ymax
    ax2.set_xlim(10**-2, 2.5)  # (0, 2.5)  # (10**-2, 13.6)
    ax2.text(2*10**-2, 4*10**-8, 'LBGs', fontsize=30)  # use if uvj_in
    # ax2.text(0.5, 5*10**-8, 'LBGs', fontsize=30)  # use if uvj_in
    '''
    if priors is not None:
        ax1.axvline(x=priors[0][1], color='k', label='Prior')  # median, +/- 1sigma for EELG prior
        ax2.axvline(x=priors[1][1], color='k', label='Prior')  # median, +/- 1sigma for LBG prior

    fs = 20
    ax1.legend(numpoints=1, loc='upper left', prop={'size': fs})
    ax2.legend(numpoints=1, loc='upper left', prop={'size': fs})
    plt.setp(ax2.get_yticklabels(), visible=False)  # hide y-axis labels on right-hand subplot to prevent overlap
    # plt.subplots_adjust(wspace=0.05)  # vertical whitespace (i.e. the width) between the two subplots
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size': 22})
    # plt.xlabel('SSFR (most recent bin)', fontsize=30)
    fig.text(0.5, 0.04, 'SSFR (most recent bin)', ha='center', fontsize=30)  # 30
    # plt.tight_layout()

    if save:
        plt.savefig(title + '.png', bbox_inches='tight')
    else:
        plt.show()


if __name__ == "__main__":

    comp = 1
    boot = 0
    vary = 1
    mask = 0
    others = 0
    short = 0
    if vary:
        base = ['vary', 'vary']
        folders = ['pkl_evar/', 'pkl_nvar/']
        mass = [10.15, 10.51]
        import eelg_varymet_params as e_params
        import eelg_varymet_params as n_params
    elif mask:
        base = ['newmask', 'newmask']
        folders = ['pkl_emask/', 'pkl_nmask/']
        mass = [10.15, 10.51]
        import eelg_newmask_params as e_params
        import eelg_newmask_params as n_params
    elif others:
        base = ['thirty', 'nth']  # base = ['otherbins', 'nother']  # use for otherbins
        folders = ['etpkls/', 'ntpkls/']  # ['opkls/', 'nopkls/']
        mass = [10., 10.3]  # PLACEHOLDER
        import eelg_thirty_params as e_params
        import noelg_thirty_params as n_params
    elif short:
        base = ['short', 'short']
        folders = ['pkl_eshort/', 'pkl_nshort/']
        mass = [10.07, 10.39]
        import eelg_short_params as e_params
        import eelg_short_params as n_params
    else:
        base = ['fixedmet', 'noelg']  # use for fixedmet
        folders = ['pkls/', 'nmpkls/']
        mass = [9.98, 10.26]
        import eelg_fixedmet_params as e_params
        import noelg_multirun_params as n_params

    pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[0]
    if comp == 0:
        eelg_list = open('eelg_specz_ids', 'r')
        eelgs = []
        e_objs = []
        e_fields = []
        for line in eelg_list:
            if line[0] == '#':
                pass
            else:
                cols = line.split()
                e_objs.append(cols[1])
                e_fields.append(cols[0])
                eelgs.append(cols[1] + '_' + cols[0] + '_' + base[0])  # base[0] = fixedmet (or otherbins)
        eelg_list.close()
    else:
        eelg_list = open('Comp_10.dat', 'r')
        eelgs = []
        for line in eelg_list:
            if line[0] == '#':
                pass
            else:
                cols = line.split()
                if int(cols[0]) - 200000 > 0:
                    eelgs.append(str(int(cols[0]) - 200000) + '_uds_' + base[0])  # base[1] = noelg (or nother)
                elif int(cols[0]) - 100000 > 0:
                    eelgs.append(str(int(cols[0]) - 100000) + '_cosmos_' + base[0])  # base[1] = noelg (or nother)
                else:
                    eelgs.append(str(int(cols[0])) + '_cdfs_' + base[0])  # base[1] = noelg (or nother)
        eelg_list.close()

    lbg_list = open('lbg_ids', 'r')
    flist = {}
    lbgs = []
    l_objs = []
    l_fields = []
    l_pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[1]
    for line in lbg_list:
        if int(line) - 200000 > 0:
            flist[str(int(line) - 200000)] = 'uds'
            lbgs.append(str(int(line) - 200000) + '_uds_' + base[1])  # base[1] = noelg (or nother)
            l_objs.append(int(line) - 200000)
            l_fields.append('uds')
        elif int(line) - 100000 > 0:
            flist[str(int(line) - 100000)] = 'cosmos'
            lbgs.append(str(int(line) - 100000) + '_cosmos_' + base[1])
            l_objs.append(int(line) - 100000)
            l_fields.append('cosmos')
        else:
            flist[str(int(line))] = 'cdfs'
            lbgs.append(str(int(line)) + '_cdfs_' + base[1])
            l_objs.append(int(line))
            l_fields.append('cdfs')
    lbg_list.close()

    tun = True
    if tun:  # flat ssfr prior = 1 / t_univ, based on perc of t_univ values for each population
        pri = [0.41772065 * 1e-9, 0.50135904 * 1e-9, 0.55399038 * 1e-9]  # USE
        pri_l = [0.40128419 * 1e-9, 0.44860297 * 1e-9, 0.56183993 * 1e-9]  # from ^ commented out thing  # USE
        # --> t_univ prior_l = [1.78*1e9, 2.23*1e9, 2.49*1e9]  # (2.23 Gyr +0.25 Gyr / -0.45 Gyr)
    else:
        pri = draw_ssfr_from_prior(['1824'], ['cosmos'], ndraw=1e4, alpha_sfh=1.0, pfile=e_params, show=False,
                                   t=1.99*1e9)  # if short, t=1e9
        pri_l = draw_ssfr_from_prior(['5957'], ['uds'], ndraw=1e4, alpha_sfh=1.0, pfile=n_params, show=False,
                                     t=2.23*1e9)  # if short, t =1.23*1e9
        print(pri)
        print(pri_l)

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
            # temp = randraw(file)  # temp[0] lists num=1000 random posterior samples; temp[1] = time vector
            temp = randraw(file, mass[0])  # temp[0] lists num=1000 random posterior samples; temp[1] = time vector
            # temp = bootdraw(file)  # temp[0] lists num=1000 random posterior samples; temp[1] = time vector
            draws.append(temp[0])
            # boots.append(bootstrap(temp[0]))
            t1.append(temp[1])
        else:
            print(file)

    # stacker2(draws, t1)
    sig = 1  # what sigma error to show on plot
    all1 = stacker(draws, sigma=sig)

    draws2 = []
    numl = 0
    cl = 0
    # t2 = []
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

    all2 = stacker(draws2, sigma=sig)

    newall1 = all1[0] + all1[1] + all1[2]
    newall2 = all2[0] + all2[1] + all2[2]
    print(len(all1[0]), len(all2[0]), len(all1[1]))

    # smooth_percs = perc1, perc2
    print(nummy, c, 'nume')
    print(numl, cl, 'numl')
    # smooth_percs = [smooth(perc1), smooth(perc2)]

    plot_sfhs([newall1, newall2], t1[0], elist=eelgs, llist=lbgs, uvj_in=True, sigma=sig, priors=[pri, pri_l], tuniv=tun)