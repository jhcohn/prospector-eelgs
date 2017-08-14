import pickle
import matplotlib.pyplot as plt
import numpy as np


def plop(file_list):
    for name in file_list:
        file = name[0] + '_' + name[1] + '_' + name[2] + '_extra_out.pkl'

        with open(file, 'rb') as data:
            extra_output = pickle.load(data)
            plt.loglog(extra_output['extras']['t_sfh'], extra_output['bfit']['sfh'], lw=2, label=name[0], alpha=0.6)
        plt.ylabel(r'Best-fit SFH [M$_\odot$ yr$^{-1}$]')
        plt.ylim(10**-3, 10**4)
        plt.xlim(10**-2, 13.6)

    plt.xlabel('Lookback time [Gyr]')
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.legend(numpoints=1, loc=0, prop={'size': 20})
    plt.rcParams.update({'font.size': 22})
    plt.show()


def stack(file_list, label):
    files = []
    for name in file_list:
        files.append(name[0] + '_' + name[1] + '_' + name[2] + '_extra_out.pkl')
    sfh = [None] * len(file_list)
    ts = [None] * len(file_list)

    for i in range(len(files)):  # same as len(file_list)
        with open(files[i], 'rb') as data:
            extra_output = pickle.load(data)
            sfh[i] = extra_output['bfit']['sfh']  # sfh is now a list where each element is an SFH array
            ts[i] = extra_output['extras']['t_sfh']

    t = ts[0]  # when bins are same, all t are the same, so doesn't matter which extra_output t_sfh we use
    avg = []
    for j in range(len(sfh[0])):
        sumsfh = 0
        for i in range(len(sfh)):
            sumsfh += sfh[i][j]  # takes sfh[0][0] + sfh[1][0] + sfh[2][0], next line avgs; repeat for sfh[0][1] + ...
        avg.append(sumsfh / len(sfh))

    plt.loglog(t, avg, lw=2, alpha=0.6)  # label='stack',
    plt.ylabel(r'Best-fit SFH [M$_\odot$ yr$^{-1}$]')
    plt.ylim(10**-3, 10**4)
    plt.xlim(10**-3, 13.6)

    plt.xlabel('Lookback time [Gyr]')
    plt.title(label)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    # plt.legend(numpoints=1, loc=0, prop={'size': 20})
    plt.rcParams.update({'font.size': 22})
    plt.show()


def dbl_stack(eelgs, lbgs):
    fig = plt.figure()
    ax1 = plt.subplot(1, 2, 1)
    ax2 = plt.subplot(1, 2, 2, sharey=ax1, sharex=ax1)

    files = []
    for name in eelgs:
        files.append(name[0] + '_' + name[1] + '_' + name[2] + '_extra_out.pkl')
    sfh = [None] * len(eelgs)
    ts = [None] * len(eelgs)

    for i in range(len(files)):  # same as len(file_list)
        with open(files[i], 'rb') as data:
            extra_output = pickle.load(data)
            sfh[i] = extra_output['bfit']['sfh']  # sfh is now a list where each element is an SFH array
            ts[i] = extra_output['extras']['t_sfh']

    t = ts[0]  # when bins are same, all t are the same, so doesn't matter which extra_output t_sfh we use
    avg = []
    for j in range(len(sfh[0])):
        sumsfh = 0
        for i in range(len(sfh)):
            sumsfh += sfh[i][j]  # takes sfh[0][0] + sfh[1][0] + sfh[2][0], next line avgs; repeat for sfh[0][1] + ...
        avg.append(sumsfh / len(sfh))

    ax1.plot(t, avg, lw=2, alpha=0.6)  # label='stack',
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_ylabel(r'Best-fit SFH [M$_\odot$ yr$^{-1}$]')
    ax1.set_ylim(10**-3, 10**4)
    ax1.set_xlim(10**-2, 13.6)

    # NOW REPEAT FOR LBGs
    files_lbg = []
    for name in lbgs:
        files_lbg.append(name[0] + '_' + name[1] + '_' + name[2] + '_extra_out.pkl')
    sfh_l = [None] * len(lbgs)
    ts_l = [None] * len(lbgs)

    for i in range(len(files_lbg)):  # same as len(file_list)
        with open(files_lbg[i], 'rb') as data:
            extra_output = pickle.load(data)
            sfh_l[i] = extra_output['bfit']['sfh']  # sfh is now a list where each element is an SFH array
            ts_l[i] = extra_output['extras']['t_sfh']

    t_l = ts_l[0]  # when bins are same, all t are the same, so doesn't matter which extra_output t_sfh we use
    avg_l = []
    for j in range(len(sfh_l[0])):
        sumsfh = 0
        for i in range(len(sfh_l)):
            sumsfh += sfh_l[i][j]  # takes sfh[0][0] + sfh[1][0] + sfh[2][0], next line avgs; repeat for sfh[0][1] + ...
        avg_l.append(sumsfh / len(sfh_l))

    ax2.plot(t_l, avg_l, lw=2, alpha=0.6)  # label='stack',
    # ax2.set_xscale('log')
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.subplots_adjust(wspace=0.05)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size': 22})
    fig.text(0.5, 0.04, 'Lookback time [Gyr]', ha='center')
    # plt.xlabel('Lookback time [Gyr]')
    # plt.legend(numpoints=1, loc=0, prop={'size': 20})
    plt.show()

if __name__ == "__main__":
    f = ['cdfs', 'cosmos', 'uds']  # ZFOURGE fields
    b = ['fixedmet', 'noelg', 'noelgduston']  # param file bases

    # EELGs
    c1824 = ['1824', f[1], b[0]]
    c12105 = ['12105', f[1], b[0]]
    c16067 = ['16067', f[1], b[0]]
    c17423 = ['17423', f[1], b[0]]
    u5206 = ['5206', f[2], b[0]]
    eelgs = [c1824, c16067, c12105, c17423, u5206]

    # STAR-FORMING NON-EELGs
    c5029 = ['5029', f[1], b[2]]
    c7730 = ['7730', f[1], b[0]]
    sfgs = [c5029, c7730]

    # QUIESCENTS
    c13110 = ['13110', f[1], b[1]]
    f20752 = ['20752', f[0], b[1]]
    quis = [c13110, f20752]

    '''
    # plop(eelgs)  # plots all the sfhs at once
    stack(eelgs, 'EELGs')  # plots avg sfr in each time bin; stack eelgs
    stack(quis, 'Quiescent galaxies')  # stack quiescent galaxies
    stack(sfgs, 'Star-forming galaxies')  # stack star-forming non-eelgs
    '''
    dbl_stack(eelgs, quis)
'''
Currently running with:
python stack_sfh.py
'''
