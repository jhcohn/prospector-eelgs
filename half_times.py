import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import random
import os
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import uvj

if __name__ == "__main__":

    others = False
    logscale = True
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

    # START STACKING
    sig = 1  # what sigma error to show on plot

    ht = []
    nummy = 0
    c = 0
    for glxy in eelgs:
        c += 1
        file = pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            nummy += 1
            # open, load extra_output
            with open(file, 'rb') as f:
                extra_output = pickle.load(f)
            # select half_time from extra_output
            ht.append(extra_output['bfit']['half_time'])
        else:
            print(file)
    print(nummy, c, 'nume')
    percs = np.percentile(ht, [16, 50, 84])

    htl = []
    numl = 0
    cl = 0
    for glxy in lbgs:
        cl += 1
        file = l_pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            numl += 1
            # open, load extra_output
            with open(file, 'rb') as f:
                extra_output = pickle.load(f)
            # select half_time from extra_output
            htl.append(extra_output['bfit']['half_time'])
        else:
            print(file)
    print(numl, cl, 'numl')
    percsl = np.percentile(htl, [16, 50, 84])

    # plot!
    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)
    if logscale:
        ax1.set_xscale('log')
    ax1.hist(ht, bins=100, weights=[1./nummy]*len(ht), histtype="step", normed=False, color='b', lw=2)
    ax1.hist(htl, bins=100, weights=[1./numl]*len(htl), histtype="step", normed=False, color='r', lw=2)

    # plot median, +/-1sigma for both histograms
    ax1.axvline(x=percs[1], color='b', linestyle='--', lw=2)
    ax1.axvline(x=percsl[1], color='r', linestyle='--', lw=2)

    # shade in +/-1sigma region
    ax1.axvspan(percs[0], percs[2], color='b', alpha=0.2)
    ax1.axvspan(percsl[0], percsl[2], color='r', alpha=0.2)

    # figure labels
    fs = 20
    hi = 0.4
    if others:
        hi = 0.25
    ax1.text(2, hi, 'EELGs', color='b', fontsize=fs)
    ax1.text(2, hi-0.02, 'LBGs', color='r', fontsize=fs)
    ax1.set_xlabel('Half-mass assembly time [Gyr]', ha='center', fontsize=fs)  # 30
    ax1.set_ylabel(r'Fraction of galaxies in sample', fontsize=fs)
    '''
    ax1 = plt.subplot(1, 2, 1)
    ax2 = plt.subplot(1, 2, 2, sharey=ax1)
    if logscale:
        ax1.set_xscale('log')
        ax2.set_xscale('log')

    # weights s.t. total in all the bins sums up to 1!
    ax1.hist(ht, bins=100, weights=[1./nummy]*len(ht), histtype="stepfilled", normed=False, lw=2)
    ax2.hist(htl, bins=100, weights=[1./numl]*len(htl), histtype="stepfilled", normed=False, lw=2)

    # plot median, +/-1sigma for both histograms
    ax1.axvline(x=percs[1], color='k', linestyle='--', lw=2)
    ax2.axvline(x=percsl[1], color='k', linestyle='--', lw=2)

    # shade in +/-1sigma region
    ax1.axvspan(percs[0], percs[2], color='k', alpha=0.2)
    ax2.axvspan(percsl[0], percsl[2], color='k', alpha=0.2)

    # figure labels
    fs = 20
    hi = 0.4
    if others:
        hi = 0.25
    ax1.text(2, hi, 'EELGs', fontsize=fs)
    ax2.text(2, hi, 'LBGs', fontsize=fs)
    fig.text(0.5, 0.04, 'Half-mass assembly time [Gyr]', ha='center', fontsize=fs)  # 30
    ax1.set_ylabel(r'Fraction of galaxies in sample', fontsize=fs)
    '''
    plt.show()
