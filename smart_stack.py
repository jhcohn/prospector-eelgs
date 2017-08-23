import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import random


def randraw(infile, num=1000):
    """
    For a given galaxy, randraw samples the posterior for each point in extra_output['extras']['sfh'][i] num times

    :param infile:
    :param num:
    :return:
    """
    with open(infile, 'rb') as exout:
        extra_output = pickle.load(exout)

    # print(len(extra_output['extras']['sfh']), len(extra_output['extras']['sfh'][0]))  # 22, 2000
    draw_from_sfh = np.zeros(shape=(len(extra_output['extras']['sfh']), num))  # shape=(22, num)

    # len(draw_from_sfh)=22; at each of these 22 points, randomly draw from the sfh posterior num=1000 times:
    for i in range(len(extra_output['extras']['sfh'])):
        for j in range(1000):
            draw_from_sfh[i][j] = extra_output['extras']['sfh'][i][random.randint(0, num)]

    return draw_from_sfh, extra_output['extras']['t_sfh']


def smooth(perc):
    """
    Takes the stacked sfh that is output from stacker and averages the sfr values in each bin, such that the sfr within
    each bin is flat, as is the case in the original extra_output['extras']['sfh'] output

    :param perc:
    :return:
    """
    # from perc: bin1 0:2, bin2 3:6, bin3 7:10, bin4 11:14, bin5 15:18, bin6 19:22
    smoother = np.zeros(shape=(len(perc), 3))
    # print(len(perc[0]))  # 22
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

    :param gal_draws:
    :return:
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
        perc[jj, :] = np.percentile(all_draws[jj, :], [16.0, 50.0, 84.0])  # median, +/- 34

    return perc  # len(perc) = 22, len(perc[0]) = 3


def plot_sfhs(percs, t, lw=1):

    fig = plt.figure()
    ax1 = plt.subplot(1, 2, 1)
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_ylim(10**-2, 10**3)
    ax1.set_xlim(10**-2, 13.6)
    ax1.set_ylabel(r'Stacked SFH [M$_\odot$ yr$^{-1}$]')
    ax1.text(4, 10**2.5, 'EELGs', fontsize=30)

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
    ax2.set_ylim(10**-2, 10**3)
    ax2.set_xlim(10**-2, 13.6)
    ax2.text(4, 10**2.5, 'LBGs', fontsize=30)

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
    eelgs = [c1824, c16067, c12105, c17423, u5206]

    # LBGs
    # c5029 = ['5029', f[1], b[2]]  # z~2.14, emission comment --> not typical necessarily
    # c7730 = ['7730', f[1], b[0]]  # z~2.1973, "gorgeous" emission --> not typical probably
    f6900 = ['6900', f[0], b[0]]
    f29430 = ['29430', f[0], b[0]]
    lbgs = [f6900, f29430]

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

    draws_lbg = []
    t_l = []
    for glxy in lbgs:
        file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
        temp = randraw(file)
        draws_lbg.append(temp[0])
        t_l.append(temp[1])

    perc_l = stacker(draws_lbg)

    percs = [smooth(perc1), smooth(perc_l)]
    plot_sfhs(percs, t[0])


'''
# run from command line in snow environment using:
python smart_stack.py
'''

'''
Vy Skype:
How to choose which LBGs to use for comparison? Just get galaxies in roughly same redshift range?
What exact redshift range to use? IDs Ben sent me ages ago were at z~2.95-3.65; want generally 3-4?
**Email Ben for IDs in composite that most resembles LBGs
'''