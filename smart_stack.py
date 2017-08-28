import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import random


def randraw(infile, num=1000):
    """
    For a given galaxy, randraw samples the posterior for each point in extra_output['extras']['sfh'][i] num times

    :param infile: ID_field_base_extra_out.py, where extra_output is stored using output.py
    :param num: number of times to sample the galaxy posterior at each point in the SFH
    :return: draw_from_sfh = 22 x num, lists the num random posterior samples, t_sfh = time vector associated with SFH
    """
    with open(infile, 'rb') as exout:
        extra_output = pickle.load(exout)

    # print(len(extra_output['extras']['ssfr']), len(extra_output['extras']['ssfr'][0]))  # 22, 2000
    draw_from_sfh = np.zeros(shape=(len(extra_output['extras']['ssfr']), num))  # shape=(22, num)
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


def stacker(gal_draws):
    """
    stacker takes input of random points drawn from a list of galaxies' SFH posteriors, concatenates them within each
    bin, and then calculates the median and 1 sigma errors in each bin

    draws should be in format draws = [draw_from_sfh1, draw_from_sfh2, ...], where each draw_from_sfh has shape=(22,num)

    :param gal_draws: list comprised of draw_from_sfh (output from randraw) for a set of galaxies
    :return: perc = stored lists of the median and +/1 1sigma SFH values calculated from the gal_draws
    """

    # len(gal_draws) = number of galaxies in stack; len(gal_draws[0]) = 22, len(gal_draws[0][0]) = 1000
    all_draws = np.zeros(shape=(len(gal_draws[0]), len(gal_draws[0][0]) * len(gal_draws)))
    for k in range(len(gal_draws)):
        note = k * 1000  # append the 1000 values in each gal_draws[k] at each of the 22 points to all_draws
        for i in range(len(gal_draws[k])):
            for j in range(len(gal_draws[k][i])):
                all_draws[i][note+j] += gal_draws[k][i][j]
    print(len(all_draws), len(all_draws[0]))

    perc = np.zeros(shape=(len(gal_draws[0]), 3))  # len(gal_draws[0]) = 22 = len(t)
    for jj in xrange(len(gal_draws[0])):
        perc[jj, :] = np.percentile(all_draws[jj, :], [16.0, 50.0, 84.0])  # median, +/- 34% = +/- 1sigma

    return perc  # len(perc) = 22, len(perc[0]) = 3


def plot_sfhs(percs, t, lw=1, spec=True):
    """
    Plots SFH stacks for two different galaxy samples side-by-side

    :param percs: list of two smoothed percs, for two different galaxy samples, each output by smooth(perc)
    :param t: time vector output by randraw
    :param lw: line width
    :param spec: if stacking specific SFR instead of plain SFR, spec=True
    :return: plot
    """
    if spec:
        ymin, ymax = 1e-12, 1e-7
    else:
        ymin, ymax = 1e-2, 1e3

    fig = plt.figure()
    ax1 = plt.subplot(1, 2, 1)
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_ylim(ymin, ymax)
    ax1.set_xlim(10**-2, 13.6)
    ax1.set_ylabel(r'Stacked sSFH [M$_\odot$ yr$^{-1}$]')
    ax1.text(4, 10**-8, 'EELGs', fontsize=30)
    # ax1.text(4, 10**2.5, 'EELGs', fontsize=30)

    ax1.plot(t, percs[0][:, 1], '-', color='k', lw=lw)  # median
    ax1.fill_between(t, percs[0][:, 0], percs[0][:, 2], color='k', alpha=0.3)  # fill region between +/- 1sigma
    ax1.plot(t, percs[0][:, 0], '-', color='k', alpha=0.3, lw=lw)  # -1sigma
    ax1.plot(t, percs[0][:, 2], '-', color='k', alpha=0.3, lw=lw)  # +1sigma

    ax2 = plt.subplot(1, 2, 2)  # , sharey=ax1, sharex=ax1)  # don't need to share axis if plot same region & never zoom
    ax2.plot(t, percs[1][:, 1], '-', color='k', lw=lw)  # median
    ax2.fill_between(t, percs[1][:, 0], percs[1][:, 2], color='k', alpha=0.3)  # fill region between +/- 1sigma
    ax2.plot(t, percs[1][:, 0], '-', color='k', alpha=0.3, lw=lw)  # -1sigma
    ax2.plot(t, percs[1][:, 2], '-', color='k', alpha=0.3, lw=lw)  # +1sigma
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_ylim(ymin, ymax)
    ax2.set_xlim(10**-2, 13.6)
    ax2.text(4, 10**-8, 'LBGs', fontsize=30)
    # ax2.text(4, 10**2.5, 'LBGs', fontsize=30)

    plt.setp(ax2.get_yticklabels(), visible=False)  # hide y-axis labels on right-hand subplot to prevent overlap
    plt.subplots_adjust(wspace=0.05)  # vertical whitespace (i.e. the width) between the two subplots
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size': 22})
    fig.text(0.5, 0.04, 'Lookback time [Gyr]', ha='center')
    plt.show()


if __name__ == "__main__":
    f = ['cdfs', 'cosmos', 'uds']  # ZFOURGE fields
    b = ['fixedmet', 'noelg', 'noelgduston', 'dust']  # param file bases

    # EELGs
    c1824 = ['1824', f[1], b[0]]
    c12105 = ['12105', f[1], b[0]]
    c16067 = ['16067', f[1], b[0]]
    c17423 = ['17423', f[1], b[0]]
    f8941 = ['8941', f[0], b[0]]
    u5206 = ['5206', f[2], b[0]]
    eelgs = [c1824, c12105, c16067, c17423, u5206]

    # LBGs
    # c5029 = ['5029', f[1], b[2]]  # z~2.14, emission comment --> not typical necessarily
    # c7730 = ['7730', f[1], b[0]]  # z~2.1973, "gorgeous" emission --> not typical probably
    c5843 = ['5843', f[1], b[3]]
    f6900 = ['6900', f[0], b[0]]
    f15921 = ['15921', f[0], b[3]]
    f29430 = ['29430', f[0], b[0]]
    lbgs = [c5843, f6900, f15921, f29430]

    # QUIESCENTS
    c13110 = ['13110', f[1], b[1]]
    f20752 = ['20752', f[0], b[1]]
    quis = [c13110, f20752]

    # START STACKING
    t = []
    draws = []
    for glxy in eelgs:
        file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
        temp = randraw(file)
        draws.append(temp[0])
        t.append(temp[1])

    perc1 = stacker(draws)

    draws2 = []
    t2 = []
    for glxy in lbgs:
        file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
        temp = randraw(file)
        draws2.append(temp[0])
        t2.append(temp[1])

    perc2 = stacker(draws2)

    smooth_percs = [smooth(perc1), smooth(perc2)]
    plot_sfhs(smooth_percs, t[0])

'''
# run from command line in snow environment using:
python smart_stack.py
'''
