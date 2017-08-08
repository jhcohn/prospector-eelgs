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


def stack(file_list):
    files = []
    for name in file_list:
        files.append(name[0] + '_' + name[1] + '_' + name[2] + '_extra_out.pkl')
    sfh = [None] * len(file_list)

    for i in range(len(files)):  # same as len(file_list)
        with open(files[i], 'rb') as data:
            extra_output = pickle.load(data)
            sfh[i] = extra_output['bfit']['sfh']  # sfh is now a list where each element is an SFH array

    t = extra_output['extras']['t_sfh']  # all t are the same, so doesn't matter which extra_output this is

    avg = []
    for j in range(len(sfh[0])):
        sumsfh = 0
        for i in range(len(sfh)):
            sumsfh += sfh[i][j]  # takes sfh[0][0] + sfh[1][0] + sfh[2][0], next line avgs; repeat for sfh[0][1] + ...
        avg.append(sumsfh / len(sfh))

    plt.loglog(t, avg, lw=2, label='stack', alpha=0.6)
    plt.ylabel(r'Best-fit SFH [M$_\odot$ yr$^{-1}$]')
    plt.ylim(10**-3, 10**4)
    plt.xlim(10**-2, 13.6)

    plt.xlabel('Lookback time [Gyr]')
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.legend(numpoints=1, loc=0, prop={'size': 20})
    plt.rcParams.update({'font.size': 22})
    plt.show()

if __name__ == "__main__":
    f = ['cdfs', 'cosmos', 'uds']  # ZFOURGE fields
    b = ['fixedmet']  # ZFOURGE bases

    c1824 = ['1824', f[1], b[0]]
    c12105 = ['12105', f[1], b[0]]
    c16067 = ['16067', f[1], b[0]]
    c17423 = ['17423', f[1], b[0]]
    eelgs = [c1824, c16067, c12105, c17423]

    plop(eelgs)  # plots all the sfhs at once

    stack(eelgs)  # plots avg sfr in each time bin

'''
Currently running with:
python stack_sfh.py
'''
