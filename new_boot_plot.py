import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import random
import os


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
    print(resample_i)
    X_resample = X[resample_i]  # take X and use indices chosen randomly above
    if X_err != None:
        X_err_resample = X_err[resample_i]
        return X_resample, X_err_resample
    else:
        return X_resample


def randraw(infile, num=1000):  # num=1000
    """
    For a given galaxy, randraw samples the posterior for each point in extra_output['extras']['sfh'][i] num times

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


def plot_sfhs(percs, t, lw=1, spec=True, sigma=1, save=False, title=None):
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
        ymin, ymax = 1e-11, 1e-7
        label = r'Stacked sSFH [yr$^{-1}$]'

    fig = plt.figure()
    ax1 = plt.subplot(1, 2, 1)
    ax2 = plt.subplot(1, 2, 2)  # , sharey=ax1, sharex=ax1)  # don't need to share axis if plot same region & never zoom

    if sigma == 1:
        ax1.plot(t, percs[0][:, 1], '-', color='k', lw=lw)  # median
        ax1.fill_between(t, percs[0][:, 0], percs[0][:, 2], color='k', alpha=0.3)  # fill region between +/- 1sigma
        ax1.plot(t, percs[0][:, 0], '-', color='k', alpha=0.3, lw=lw)  # -1sigma
        ax1.plot(t, percs[0][:, 2], '-', color='k', alpha=0.3, lw=lw)  # +1sigma

        ax2.plot(t, percs[1][:, 1], '-', color='k', lw=lw)  # median
        ax2.fill_between(t, percs[1][:, 0], percs[1][:, 2], color='k', alpha=0.3)  # fill region between +/- 1sigma
        ax2.plot(t, percs[1][:, 0], '-', color='k', alpha=0.3, lw=lw)  # -1sigma
        ax2.plot(t, percs[1][:, 2], '-', color='k', alpha=0.3, lw=lw)  # +1sigma

    elif sigma == 3:
        if spec:
            ymin, ymax = 1e-13, 1e-6
        else:
            ymin, ymax = 1e-3, 1e4
        ax1.plot(t, percs[0][:, 3], '-', color='k', lw=lw)  # median
        ax1.fill_between(t, percs[0][:, 2], percs[0][:, 4], color='k', alpha=0.3)  # fill region between +/- 1sigma
        ax1.fill_between(t, percs[0][:, 1], percs[0][:, 5], color='k', alpha=0.3)  # fill region between +/- 2sigma
        ax1.fill_between(t, percs[0][:, 0], percs[0][:, 6], color='k', alpha=0.3)  # fill region between +/- 3sigma
        ax1.plot(t, percs[0][:, 0], '-', color='k', alpha=0.3, lw=lw)  # -3sigma
        ax1.plot(t, percs[0][:, 6], '-', color='k', alpha=0.3, lw=lw)  # +3sigma

        ax2.plot(t, percs[1][:, 3], '-', color='k', lw=lw)  # median
        ax2.fill_between(t, percs[1][:, 2], percs[1][:, 4], color='k', alpha=0.3)  # fill region between +/- 1sigma
        ax2.fill_between(t, percs[1][:, 1], percs[1][:, 5], color='k', alpha=0.3)  # fill region between +/- 2sigma
        ax2.fill_between(t, percs[1][:, 0], percs[1][:, 6], color='k', alpha=0.3)  # fill region between +/- 3sigma
        ax2.plot(t, percs[1][:, 0], '-', color='k', alpha=0.3, lw=lw)  # -3sigma
        ax2.plot(t, percs[1][:, 6], '-', color='k', alpha=0.3, lw=lw)  # +3sigma

    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_ylim(ymin, ymax)
    ax1.set_xlim(10**-2, 2.5)  # (0, 2.5)  # (10**-2, 13.6)
    ax1.set_ylabel(label)
    # ax1.text(4, 10**2.5, 'EELGs', fontsize=30)
    ax1.text(1, 4*10**-8, 'EELGs', fontsize=30)

    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_ylim(ymin, ymax)
    ax2.set_xlim(10**-2, 2.5)  # (0, 2.5)  # (10**-2, 13.6)
    # ax2.text(4, 10**2.5, 'LBGs', fontsize=30)
    ax2.text(1, 4*10**-8, 'LBGs', fontsize=30)

    plt.setp(ax2.get_yticklabels(), visible=False)  # hide y-axis labels on right-hand subplot to prevent overlap
    plt.subplots_adjust(wspace=0.05)  # vertical whitespace (i.e. the width) between the two subplots
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size': 22})
    fig.text(0.5, 0.04, 'Lookback time [Gyr]', ha='center')
    plt.tight_layout()
    if save:
        plt.savefig(title + '.png', bbox_inches='tight')
    else:
        plt.show()


if __name__ == "__main__":

    others = False
    if others:
        base = ['thirty', 'nth']  # base = ['otherbins', 'nother']  # use for otherbins
        folders = ['etpkls/', 'ntpkls/']  # ['opkls/', 'nopkls/']
    else:
        base = ['fixedmet', 'noelg']  # use for fixedmet
        folders = ['pkls/', 'nmpkls/']

    eelg_list = open('eelg_specz_ids', 'r')
    pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[0]
    eelgs = []
    for line in eelg_list:
        if line[0] == '#':
            pass
        else:
            cols = line.split()
            eelgs.append(cols[1] + '_' + cols[0] + '_' + base[0])  # base[0] = fixedmet (or otherbins)
    eelg_list.close()

    lbg_list = open('lbg_ids', 'r')
    flist = {}
    lbgs = []
    l_pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[1]
    for line in lbg_list:
        if int(line) - 200000 > 0:
            flist[str(int(line) - 200000)] = 'uds'
            lbgs.append(str(int(line) - 200000) + '_uds_' + base[1])  # base[1] = noelg (or nother)
        elif int(line) - 100000 > 0:
            flist[str(int(line) - 100000)] = 'cosmos'
            lbgs.append(str(int(line) - 100000) + '_cosmos_' + base[1])
        else:
            flist[str(int(line))] = 'cdfs'
            lbgs.append(str(int(line)) + '_cdfs_' + base[1])
    lbg_list.close()

    cut = 'max2'  # 'thirty_30'
    rangey = 4
    eelg_percs = []
    lbg_percs = []
    for num in range(rangey):
        print('loop:', num)
        eelgs = bootstrap(np.asarray(eelgs))
        lbgs = bootstrap(np.asarray(lbgs))
        eelgs = eelgs[:]  # cut
        lbgs = lbgs[:]  # cut

        t1 = []
        draws = []
        boots = []
        nummy = 0
        c = 0
        for glxy in eelgs:
            c += 1
            file = pkls + glxy + '_extra_out.pkl'
            if os.path.exists(file):
                nummy += 1
                temp = randraw(file)  # temp[0] lists the num=1000 random posterior samples; temp[1] = time vector
                draws.append(temp[0])
                t1.append(temp[1])
            else:
                print(file)

        # stacker2(draws, t1)
        sig = 1  # what sigma error to show on plot
        perc1 = stacker(draws, sigma=sig)

        draws2 = []
        numl = 0
        cl = 0
        # t2 = []
        for glxy in lbgs:
            cl += 1
            file = l_pkls + glxy + '_extra_out.pkl'
            if os.path.exists(file):
                numl += 1
                temp = randraw(file)
                draws2.append(temp[0])
            else:
                print(file)

        print(nummy, c, 'nume')
        print(numl, cl, 'numl')

        perc2 = stacker(draws2, sigma=sig)

        smooth_percs = [smooth(perc1), smooth(perc2)]  # smooth_percs[0][:, 1] (or smooth_percs[1][:, 1]) = median

        eelg_percs.append(smooth_percs[0])
        lbg_percs.append(smooth_percs[1])

        print(str(num))
        # plot_sfhs(smooth_percs, t1[0], sigma=sig, save=True, title='booties/bootstrap-' + str(cut) + '_' + str(num))

    etemp_array = np.zeros(shape=(len(smooth_percs[0]), rangey))  # num of bootstraps = rangey
    ltemp_array = np.zeros(shape=(len(smooth_percs[1]), rangey))  # num of bootstraps = rangey

    print(len(eelg_percs), len(smooth_percs[0]), len(eelg_percs[0][0]))  # (2, 22, 3)
    # FIX THINGS STARTING FROM HERE

    new_draws = []
    for perc in eelg_percs:
        new_draws.append(perc)
    newstack = stacker(new_draws)

    lnew_draws = []
    for perc in lbg_percs:
        lnew_draws.append(perc)
    lnewstack = stacker(lnew_draws)

    boot_percs = [smooth(newstack), smooth(lnewstack)]

    '''
    for i in range(len(eelg_percs)):  # for i in rangey (i.e. for each galaxy)
        for j in range(len(eelg_percs[i])):  # = len(perc) = 22 (at each bin sample point)
            for k in range(len(eelg_percs[i][j])):
                etemp_array[i][j] = eelg_percs[i][j, 1]

    for i in range(len(lbg_percs)):
        for j in range(len(lbg_percs[i])):  # = len(perc) = 22
            ltemp_array[i][j] = lbg_percs[i][j, 1]

    eboot_perc = np.zeros(shape=(len(smooth_percs[0]), 3))  # len(gal_draws[0])=22=len(t); len(perc)=22, len(perc[0])=3
    lboot_perc = np.zeros(shape=(len(smooth_percs[0]), 3))  # len(gal_draws[0])=22=len(t); len(perc)=22, len(perc[0])=3
    for jj in xrange(len(smooth_percs[0])):
        eboot_perc[jj, :] = np.percentile(etemp_array[jj, :], [16.0, 50.0, 84.0])  # median, +/- 34% = +/- 1sigma
    for jj in xrange(len(smooth_percs[0])):
        lboot_perc[jj, :] = np.percentile(etemp_array[jj, :], [16.0, 50.0, 84.0])  # median, +/- 34% = +/- 1sigma

    boot_percs = [eboot_perc, lboot_perc]
    '''
    plot_sfhs(boot_percs, t1[0], sigma=1, save=False)  # , title='booties/bootstrap-' + str(cut) + '_' + str(num))

'''
Now create new error plot with errors determined by difference in the median due to bootstrapping
Note: len(perc) = 22, len(perc[0]) = 3 (smooth_percs[0] = perc1, smoothed)
temp_array_0 = np.zeros(shape=(len(smooth_percs[0]), rangey))  # num of bootstraps = rangey
temp_array_1 = np.zeros(shape=(len(smooth_percs[1]), rangey))  # num of bootstraps = rangey
for i in range(rangey):
    for j in range(len(smooth_percs[i])):  # = len(perc) = 22
        temp_array[i, j] = smooth_percs[i][j, 1]

for jj in xrange(len(smooth_percs[0])):
    boot_perc[jj, :] = np.percentile(temp_array[jj, :], [16.0, 50.0, 84.0])  # median, +/- 34% = +/- 1sigma


'''