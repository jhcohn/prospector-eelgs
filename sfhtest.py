import prospect.io.read_results as bread
import matplotlib.pyplot as plt
import numpy as np
from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log
import sys
import argparse
import pickle
import os
import get_mass_dust as gmd
import ssfr_fig3 as sfh
from matplotlib import rc
from scipy import stats
np.errstate(invalid='ignore')


if __name__ == "__main__":

    folder = 'out_newsfhtest2/'
    # pkl = 'pkl_sfhtest/'
    # pars = ['eelg_test' + str(n) + '_params.py' for n in range(100)]

    gals = {}
    git = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
    oute = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + folder
    # pkls = git + pkl
    allkeys = []
    for file in os.listdir(oute):
        c = 0
        key = ''
        for char in file:
            if char == '_':
                c += 1
            elif c == 3:
                key += char
        if file.endswith(".h5") and float(key) < 200:
            gals[key] = file
            allkeys.append(key)

    print(len(gals), 'lens')
    print(gals)
    fig = plt.figure()
    '''
    ax1 = plt.subplot2grid(shape=(2, 6), loc=(0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)
    ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)
    ax4 = plt.subplot2grid((2, 6), (1, 1), colspan=2)
    ax5 = plt.subplot2grid((2, 6), (1, 3), colspan=2)
    '''
    ax1 = plt.subplot(2, 3, 1, aspect='equal')
    ax2 = plt.subplot(2, 3, 2, aspect='equal')
    ax3 = plt.subplot(2, 3, 3, aspect='equal')
    ax4 = plt.subplot(2, 3, 4, aspect='equal')
    ax5 = plt.subplot(2, 3, 5, aspect='equal')
    ax6 = plt.subplot(2, 3, 6, aspect='equal')

    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    fs_text = 25  # 30
    fs = 15  # 20
    fs_ticks = 20  # 25

    l1, l2, l3, l4, l5, l6 = [[], []], [[], []], [[], []], [[], []], [[], []], [[], []]

    # get_e = np.zeros(shape=(4, len(gals)))  # 6 rows (mass, dust, logzsol, gaslogz, ssfr1, ssfr2)
    #onedraw = np.zeros(shape=(80, 6, 10**3))  # *10**3))  # [mass, dust, metal, gasmet, ssfr1, ssfr2]
    onedraw = np.zeros(shape=(100, 6, 10**3))  # *10**3))  # [mass, dust, metal, gasmet, ssfr1, ssfr2]  # BUCKET
    # (shape=(len(gals), 6, 10**3))
    # NOTE: INPUT SAME gasmet each time
    # sfh_bins = np.zeros(shape=(len(gals), 2, 10**4))  # 10**4 = num arg in ssfr_fig3.randraw()
    # for ind in range(len(allkeys)):
    #offsets = np.zeros(shape=(80, 6))  # [mass, dust, metal, sfh1, sfh2, sfh1/sfh2]
    #meds = np.zeros(shape=(80, 6))  # [mass, dust, metal, sfh1, sfh2, sfh1/sfh2]
    offsets = np.zeros(shape=(100, 6))  # [mass, dust, metal, sfh1, sfh2, sfh1/sfh2]  # BUCKET
    meds = np.zeros(shape=(100, 6))  # [mass, dust, metal, sfh1, sfh2, sfh1/sfh2]  # BUCKET
    #for i in range(20, 100):  # for each galaxy (0 through 99)
    for i in range(100):  # for each galaxy (0 through 99)
        # i = allkeys[ind]
        print(i, gals[str(i)])  # i = key
        #j = i - 20  # BUCKET
        if os.path.exists(oute + gals[str(i)]):
            #onedraw[j, :, :] = gmd.printer(oute + gals[str(i)], percs=False, sfhtest=True, draw1=True)
            onedraw[i, :, :] = gmd.printer(oute + gals[str(i)], percs=False, sfhtest=True, draw1=True)  # BUCKET

        # pars = git + 'eetest/sfhtest_' + str(i) + '_params.py'
        pars = git + 'new/sfhtest_' + str(i) + '_params.py'
        #i = j  # BUCKET
        with open(pars, 'r') as parfile:
            for line in parfile:
                if line.startswith("# Codename: "):
                    spaces = 0
                    inmass, indust, inmet, insfr1, insfr2 = '', '', '', '', ''
                    for char in line:
                        if char == ' ' or char == '_':
                            spaces += 1
                        elif spaces == 2:
                            inmass += char
                        elif spaces == 3:
                            indust += char
                        elif spaces == 4:
                            inmet += char
                        elif spaces == 5:
                            insfr1 += char
                        elif spaces == 6:
                            insfr2 += char
                elif line.startswith("# Identity: "):
                    ispaces = 0
                    zred = ''
                    for char in line:
                        if char == ' ' or char == '_':
                            ispaces += 1
                        elif ispaces == 4:
                            zred += char
        zred = float(zred)
        print(zred, 'zred')

        # print(onedraw[i, 0, :])

        # ssfr050 = np.percentile(recent, [16., 50., 84.])
        # ssfr50100 = np.percentile(second, [16., 50., 84.])
        mass = np.percentile(onedraw[i, 0, :], [16., 50., 84.])
        # print(mass, 'mass')
        dust = np.percentile(onedraw[i, 1, :], [16., 50., 84.])
        met = np.percentile(onedraw[i, 2, :], [16., 50., 84.])
        gasmet = np.percentile(onedraw[i, 3, :], [16., 50., 84.])
        sfr1 = np.percentile(onedraw[i, 4, :], [16., 50., 84.])
        sfr2 = np.percentile(onedraw[i, 5, :], [16., 50., 84.])
        ratio = []
        for k in range(10**3):
            ratio.append(onedraw[i, 4, np.random.randint(0, 999)] / onedraw[i, 5, np.random.randint(0, 999)])
        sfr_ratio = np.percentile(ratio, [16., 50., 84.])
        offsets[i, 0] = mass[1] - float(inmass)
        offsets[i, 1] = dust[1] - float(indust)
        offsets[i, 2] = 10**met[1] - 10**float(inmet)
        offsets[i, 3] = sfr1[1] - float(insfr1)
        offsets[i, 4] = sfr2[1] - float(insfr2)
        offsets[i, 5] = sfr1[1] / sfr2[1] - float(insfr1) / float(insfr2)

        meds[i, 0] = mass[1]
        meds[i, 1] = dust[1]
        meds[i, 2] = 10**met[1]
        meds[i, 3] = sfr1[1]
        meds[i, 4] = sfr2[1]
        meds[i, 5] = sfr1[1] / sfr2[1]

        color = 'purple'

        ax1.errorbar(float(inmass), mass[1], yerr=np.array([[mass[1] - mass[0], mass[2] - mass[1]]]).T, fmt='o',
                     color=color)
        ax2.errorbar(float(indust) / 1.086, dust[1] / 1.086,
                     yerr=np.array([[(dust[1] - dust[0]) / 1.086, (dust[2] - dust[1]) / 1.086]]).T, fmt='o',
                     color=color)
        ax3.errorbar(10**float(inmet), 10**met[1], yerr=np.array([[10**met[1]-10**met[0], 10**met[2]-10**met[1]]]).T,
                     fmt='o', color=color)
        ax4.errorbar(float(insfr1), sfr1[1], yerr=np.array([[sfr1[1] - sfr1[0], sfr1[2] - sfr1[1]]]).T, fmt='o',
                     color=color)
        ax5.errorbar(float(insfr2), sfr2[1], yerr=np.array([[sfr2[1] - sfr2[0], sfr2[2] - sfr2[1]]]).T, fmt='o',
                     color=color)
        ax6.errorbar(float(insfr1) / float(insfr2), sfr_ratio[1],
                     yerr=np.array([[sfr_ratio[1] - sfr_ratio[0], sfr_ratio[2] - sfr_ratio[1]]]).T, fmt='o',
                     color=color)
        # '''
        l1[0].append(float(inmass))
        l1[1].append(mass[1])
        l2[0].append(float(indust))
        l2[1].append(dust[1])
        l3[0].append(10**float(inmet))
        l3[1].append(10**met[1])
        l4[0].append(float(insfr1))
        l4[1].append(sfr1[1])
        l5[0].append(float(insfr2))
        l5[1].append(sfr2[1])
        l6[0].append(float(insfr1) / float(insfr2))
        l6[1].append(sfr1[1] / sfr2[1])
        print(float(inmass), mass[1], pars)
        print(float(insfr1), sfr1[1])
        print(float(insfr2), sfr2[1])

        # print('ssfr 0-50', ssfr050[1], '+', ssfr050[2] - ssfr050[1], '-', ssfr050[1]-ssfr050[0])
        # print('ssfr 50-100', ssfr50100[1], '+', ssfr50100[2] - ssfr50100[1], '-', ssfr50100[1]-ssfr50100[0])
        print('mass', mass[1], '+', mass[2]-mass[1], '-', mass[1]-mass[0])
        print('dust', dust[1] / 1.086, '+', (dust[2]-dust[1]) / 1.086, '-', (dust[1]-dust[0]) / 1.086)
        print('met', 10**met[1], '+', 10**met[2]-10**met[1], '-', 10**met[1]-10**met[0])
        print('gasmet', 10**gasmet[1], '+', 10**gasmet[2]-10**gasmet[1], '-', 10**gasmet[1]-10**gasmet[0])

    '''
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs, pad=10)
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    plt.tick_params(axis='y', which='minor')
    ax1.set_xticks([10 ** 3, 2 * 10 ** 3, 5 * 10 ** 3, 10 ** 4, 2 * 10 ** 4])  # technically works
    ax1.set_xticklabels([r'$10^3$', r'$2\times10^3$', r'$5 \times 10^3$', r'$10^4$', r'$2\times10^4$'], size=fs_ticks)
    '''

    '''
    # EELG-LIKE RECOVERY
    efold = 'out_newsfhtest/'  # 'out_sfhhometest1/'# 'out_eetest2/'  # out_eetest
    egals = {}
    # eoute = '/home/jonathan/' + efold # '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + efold
    eoute = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + efold
    # pkls = git + pkl
    ekeys = []
    for file in os.listdir(eoute):
        c = 0
        key = ''
        for char in file:
            if char == '_':
                c += 1
            elif c == 3:
                key += char
        if file.endswith(".h5") and float(key) < 200:
            egals[key] = file
            ekeys.append(key)
    # NOTE: for eetest2 these three all had len(12, 6[, 10**3])  # len(1, 6[, 10**3])
    eonedraw = np.zeros(shape=(1, 6, 10 ** 3))  # *10**3))  # [mass, dust, metal, gasmet, ssfr1, ssfr2]
    eoffsets = np.zeros(shape=(1, 6))  # [mass, dust, metal, sfh1, sfh2, sfh1/sfh2]
    emeds = np.zeros(shape=(1, 6))  # [mass, dust, metal, sfh1, sfh2, sfh1/sfh2]
    for i in range(1):# (len(egals)):  # (12): # for each galaxy (0 through 99)
        # i = allkeys[ind]
        print(i, egals[str(i)])  # i = key
        if os.path.exists(eoute + egals[str(i)]):
            eonedraw[i, :, :] = gmd.printer(eoute + egals[str(i)], percs=False, sfhtest=True, draw1=True)  # BUCKET BREAK
            # eonedraw[i, :, :] = gmd.printer(eoute + egals[str(i)], percs=False, sfhtest=True, draw1=True)
        # pars = git + 'eetest2/sfhtest_' + str(i) + '_params.py'  # eetest2
        pars = git + 'new/sfhtest_' + str(i) + '_params.py'  # 'eehometest/sfhtest_' + str(i) + '_params.py'  # eetest2
        with open(pars, 'r') as parfile:
            for line in parfile:
                if line.startswith("# Codename: "):
                    spaces = 0
                    inmass, indust, inmet, insfr1, insfr2 = '', '', '', '', ''
                    for char in line:
                        if char == ' ' or char == '_':
                            spaces += 1
                        elif spaces == 2:
                            inmass += char
                        elif spaces == 3:
                            indust += char
                        elif spaces == 4:
                            inmet += char
                        elif spaces == 5:
                            insfr1 += char
                        elif spaces == 6:
                            insfr2 += char
                elif line.startswith("# Identity: "):
                    ispaces = 0
                    zred = ''
                    for char in line:
                        if char == ' ' or char == '_':
                            ispaces += 1
                        elif ispaces == 4:
                            zred += char
        zred = float(zred)
        print(zred, 'zred')
        mass = np.percentile(eonedraw[i, 0, :], [16., 50., 84.])
        # print(mass, 'mass')
        dust = np.percentile(eonedraw[i, 1, :], [16., 50., 84.])
        met = np.percentile(eonedraw[i, 2, :], [16., 50., 84.])
        gasmet = np.percentile(eonedraw[i, 3, :], [16., 50., 84.])
        sfr1 = np.percentile(eonedraw[i, 4, :], [16., 50., 84.])
        sfr2 = np.percentile(eonedraw[i, 5, :], [16., 50., 84.])
        ratio = []
        for k in range(10**3):
            ratio.append(eonedraw[i, 4, np.random.randint(0, 999)] / eonedraw[i, 5, np.random.randint(0, 999)])
        sfr_ratio = np.percentile(ratio, [16., 50., 84.])
        eoffsets[i, 0] = mass[1] - float(inmass)
        eoffsets[i, 1] = dust[1] - float(indust)
        eoffsets[i, 2] = 10 ** met[1] - 10 ** float(inmet)
        eoffsets[i, 3] = sfr1[1] - float(insfr1)
        eoffsets[i, 4] = sfr2[1] - float(insfr2)
        eoffsets[i, 5] = sfr1[1] / sfr2[1] - float(insfr1) / float(insfr2)
        emeds[i, 0] = mass[1]
        emeds[i, 1] = dust[1]
        emeds[i, 2] = 10 ** met[1]
        emeds[i, 3] = sfr1[1]
        emeds[i, 4] = sfr2[1]
        emeds[i, 5] = sfr1[1] / sfr2[1]
        color = 'purple'
        fmt = '*'
        sz = 20
        ax1.errorbar(float(inmass), mass[1], yerr=np.array([[mass[1] - mass[0], mass[2] - mass[1]]]).T, fmt='*',
                     color=color, markersize=sz)
        ax2.errorbar(float(indust) / 1.086, dust[1] / 1.086,
                     yerr=np.array([[(dust[1] - dust[0])/1.086, (dust[2] - dust[1])/1.086]]).T, fmt='*',
                     color=color, markersize=sz)
        ax3.errorbar(10 ** float(inmet), 10 ** met[1],
                     yerr=np.array([[10 ** met[1] - 10 ** met[0], 10 ** met[2] - 10 ** met[1]]]).T,
                     fmt='*', color=color, markersize=sz)
        ax4.errorbar(float(insfr1), sfr1[1], yerr=np.array([[sfr1[1] - sfr1[0], sfr1[2] - sfr1[1]]]).T, fmt='*',
                     color=color, markersize=sz)
        ax5.errorbar(float(insfr2), sfr2[1], yerr=np.array([[sfr2[1] - sfr2[0], sfr2[2] - sfr2[1]]]).T, fmt='*',
                     color=color, markersize=sz)
        ax6.errorbar(float(insfr1) / float(insfr2), sfr_ratio[1],
                     yerr=np.array([[sfr_ratio[1] - sfr_ratio[0], sfr_ratio[2] - sfr_ratio[1]]]).T, fmt='*',
                     color=color, markersize=sz)
        '''

    # SET UP AXES:
    axes = [ax1, ax2, ax3, ax4, ax5]
    for ax in axes:
        ax.tick_params('x', length=3, width=1, which='both', labelsize=fs, pad=10)
        ax.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)

    things = ['mass', 'dust', 'met', 'ssf1','ssfr2', 'ratio']
    offs = []
    for j in range(len(offsets[0])):  # for all 5 recorded parameter offsets
        print(np.percentile(offsets[:, j], [16., 50., 84]), 'percentiles offset')
        print(np.mean(offsets[:, j]), things[j], 'mean offset')
        offs.append(np.mean(offsets[:, j]))

    m, b = np.polyfit(l1[0], l1[1], 1)
    print(m, b, 'mass')
    diff = []
    for j in range(len(l1[0])):
        expected = l1[0][j] * m + b
        diff.append(l1[1][j] - expected)
    perc = np.percentile(diff, [16., 50., 84.])
    print('mean offset = ' + str(offs[0]), 'scatter = ' + str(perc[2] - perc[0]))

    m, b = np.polyfit(l2[0], l2[1], 1)
    print(m, b, 'dust')
    diff = []
    for j in range(len(l2[0])):
        expected = l2[0][j] * m + b
        diff.append(l2[1][j] - expected)
    perc = np.percentile(diff, [16., 50., 84.])
    print('mean offset = ' + str(offs[1]), 'scatter = ' + str(perc[2] - perc[0]))

    m, b = np.polyfit(l3[0], l3[1], 1)
    print(m, b, 'met')
    diff = []
    for j in range(len(l3[0])):
        expected = l3[0][j] * m + b
        diff.append(l3[1][j] - expected)
    perc = np.percentile(diff, [16., 50., 84.])
    print('mean offset = ' + str(offs[2]), 'scatter = ' + str(perc[2] - perc[0]))

    m, b = np.polyfit(l4[0], l4[1], 1)
    print(m, b, 'sfr1')
    diff = []
    for j in range(len(l4[0])):
        expected = l4[0][j] * m + b
        diff.append(l4[1][j] - expected)
    perc = np.percentile(diff, [16., 50., 84.])
    print('mean offset = ' + str(offs[3]), 'scatter = ' + str(perc[2] - perc[0]))

    m, b = np.polyfit(l5[0], l5[1], 1)
    print(m, b, 'sfr2')
    diff = []
    for j in range(len(l5[0])):
        expected = l5[0][j] * m + b
        diff.append(l5[1][j] - expected)
    perc = np.percentile(diff, [16., 50., 84.])
    print('mean offset = ' + str(offs[4]), 'scatter = ' + str(perc[2] - perc[0]))

    m, b = np.polyfit(l6[0], l6[1], 1)
    print(m, b, 'sfr ratio')
    diff = []
    for j in range(len(l6[0])):
        expected = l6[0][j] * m + b
        diff.append(l6[1][j] - expected)
    perc = np.percentile(diff, [16., 50., 84.])
    print('mean offset = ' + str(offs[5]), 'scatter = ' + str(perc[2] - perc[0]))

    # MASS AXES
    ax1.set_xlim(8., 10.5)
    ax1.set_ylim(8., 10.5)
    # ax1.set_xlim(8., 11.)
    # ax1.set_ylim(8., 11.)
    ax1.plot([8., 11.], [8., 11.], color='k', linestyle='--')
    ax1.set_xlabel(r'Input $\log_{10}($M/M$_\odot$)', fontsize=fs_text)
    ax1.set_ylabel(r'Recovered $\log_{10}($M/M$_\odot$)', fontsize=fs_text)
    ax1.set_xticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5])  # technically works
    ax1.set_xticklabels([r'$8.0$', r'$8.5$', r'$9.0$', r'$9.5$', r'$10.0$', r'$10.5$'], size=fs_ticks)
    ax1.set_yticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5])  # technically works
    ax1.set_yticklabels([r'$8.0$', r'$8.5$', r'$9.0$', r'$9.5$', r'$10.0$', r'$10.5$'], size=fs_ticks)
    '''
    ax1.set_xticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0])  # technically works
    ax1.set_xticklabels([r'$8.0$', r'$8.5$', r'$9.0$', r'$9.5$', r'$10.0$', r'$10.5$', r'$11.0$'], size=fs_ticks)
    ax1.set_yticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0])  # technically works
    ax1.set_yticklabels([r'$8.0$', r'$8.5$', r'$9.0$', r'$9.5$', r'$10.0$', r'$10.5$', r'$11.0$'], size=fs_ticks)
    '''

    # DUST AXES
    ax2.set_xlim(0., 1.)
    ax2.set_ylim(0., 1.)
    # ax2.set_xlim(0., 2.)
    # ax2.set_ylim(0., 2.)
    ax2.plot([0., 2.], [0., 2.], color='k', linestyle='--')
    ax2.set_xlabel(r'Input A$_V$', fontsize=fs_text)
    ax2.set_ylabel(r'Recovered A$_V$', fontsize=fs_text)
    ax2.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax2.set_xticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax2.set_yticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    '''
    ax2.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])  # technically works
    ax2.set_xticklabels([r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$'], size=fs_ticks)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])  # technically works
    ax2.set_yticklabels([r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$'], size=fs_ticks)
    '''

    # METALLICITY AXES
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlim(10 ** -2., 2*10 ** -1.)
    ax3.set_ylim(10 ** -2., 2*10 ** -1.)
    # ax3.set_xlim(10**-2., 1.)
    # ax3.set_ylim(10**-2., 1.)
    ax3.plot([10**-2., 1.], [10**-2., 1.], color='k', linestyle='--')
    # ax3.set_xlabel(r'Input $\log_{10}($Z$^*/$Z$_\odot$)', fontsize=fs_text)
    # ax3.set_ylabel(r'Recovered $\log_{10}($Z$^*/$Z$_\odot$)', fontsize=fs_text)
    ax3.set_xlabel(r'Input Z$^*/$Z$_\odot$', fontsize=fs_text)
    ax3.set_ylabel(r'Recovered Z$^*/$Z$_\odot$', fontsize=fs_text)
    ax3.set_xticks([10 ** -2, 2*10**-2, 5*10**-2, 10 ** -1, 2*10**-1])  # technically works
    ax3.set_xticklabels([r'$0.01$', r'$0.02$', r'$0.05$', r'$0.1$', r'$0.2$'], size=fs_ticks)
    ax3.set_yticks([10 ** -2, 2*10**-2, 5*10**-2, 10 ** -1, 2*10**-1])  # technically works
    ax3.set_yticklabels([r'$0.01$', r'$0.02$', r'$0.05$', r'$0.1$', r'$0.2$'], size=fs_ticks)
    '''
    ax3.set_xticks([10 ** -2, 10 ** -1, 10 ** 0])  # technically works
    ax3.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$'], size=fs_ticks)
    ax3.set_yticks([10 ** -2, 10 ** -1, 10 ** 0])  # technically works
    ax3.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$'], size=fs_ticks)

    ax3.set_xticks([10 ** -2, 2*10**-2, 5*10**-2, 10 ** -1])  # technically works
    ax3.set_xticklabels([r'$10^{-2}$', r'$2\times10^{-2}$', r'$5\times10^{-2}$', r'$10^{-1}$'], size=fs_ticks)
    '''

    # SFR BIN1 AXES
    ax4.set_xlim(0.1, 0.7)
    ax4.set_ylim(0.1, 0.7)
    # ax4.set_xlim(0., 1.)
    # ax4.set_ylim(0., 1.)
    ax4.plot([0., 1.], [0., 1.], color='k', linestyle='--')
    ax4.set_xlabel(r'Input SFH [$0 - 50$ Myr]', fontsize=fs_text)
    ax4.set_ylabel(r'Recovered SFH [$0 - 50$ Myr]', fontsize=fs_text)
    ax4.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])  # technically works
    ax4.set_xticklabels([r'$0.1$', r'$0.2$', r'$0.3$', r'$0.4$', r'$0.5$', r'$0.6$', r'$0.7$'], size=fs_ticks)
    ax4.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])  # technically works
    ax4.set_yticklabels([r'$0.1$', r'$0.2$', r'$0.3$', r'$0.4$', r'$0.5$', r'$0.6$', r'$0.7$'], size=fs_ticks)
    '''
    ax4.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax4.set_xticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    ax4.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax4.set_yticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    '''

    # SFR BIN2 AXES
    ax5.set_xlim(0., 0.3)
    ax5.set_ylim(0., 0.3)
    # ax5.set_xlim(0., 1.)
    # ax5.set_ylim(0., 1.)
    ax5.plot([0., 1.], [0., 1.], color='k', linestyle='--')
    ax5.set_xlabel(r'Input SFH [$50 - 100$ Myr]', fontsize=fs_text)
    ax5.set_ylabel(r'Recovered SFH [$50 - 100$ Myr]', fontsize=fs_text)
    ax5.set_xticks([0.0, 0.1, 0.2, 0.3])  # technically works
    ax5.set_xticklabels([r'$0.0$', r'$0.1$', r'$0.2$', r'$0.3$'], size=fs_ticks)
    ax5.set_yticks([0.0, 0.1, 0.2, 0.3])  # technically works
    ax5.set_yticklabels([r'$0.0$', r'$0.1$', r'$0.2$', r'$0.3$'], size=fs_ticks)
    '''
    ax5.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax5.set_xticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    ax5.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax5.set_yticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    '''

    # SFR RATIO AXES
    ax6.set_xscale("log")
    ax6.set_yscale("log")
    ax6.set_xlim(0.5, 15.)
    ax6.set_ylim(0.5, 15.)
    # ax6.set_xlim(10**-2., 10**2.)
    # ax6.set_ylim(10**-2., 10**2.)
    ax6.plot([10**-2., 10**2.], [10**-2., 10**2.], color='k', linestyle='--')
    ax6.set_xlabel(r'Input SFH ratio [$f_{0-50} / f_{50-100}$]', fontsize=fs_text)
    ax6.set_ylabel(r'Recovered SFH ratio [$f_{0-50} / f_{50-100}$]', fontsize=fs_text)
    ax6.set_xticks([1., 2., 5., 10.])  # technically works
    ax6.set_xticklabels([r'$1$', r'$2$', r'$5$', r'$10$'], size=fs_ticks)
    ax6.set_yticks([1., 2., 5., 10.])  # technically works
    ax6.set_yticklabels([r'$1$', r'$2$', r'$5$', r'$10$'], size=fs_ticks)
    # ax6.set_xticks([10**0., 10**1.])  # technically works
    # ax6.set_xticklabels([r'$10^{0}$', r'$10^{1}$'], size=fs_ticks)
    # ax6.set_yticks([10**0., 10**1.])  # technically works
    # ax6.set_yticklabels([r'$10^{0}$', r'$10^{1}$'], size=fs_ticks)
    '''
    ax6.set_xticks([10**-2., 10**-1., 10**0., 10**1., 10**2.])  # technically works
    ax6.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$'], size=fs_ticks)
    ax6.set_yticks([10**-2., 10**-1., 10**0., 10**1., 10**2.])  # technically works
    ax6.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$'], size=fs_ticks)
    # ax6.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    # ax6.set_yticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    '''
    plt.show()

