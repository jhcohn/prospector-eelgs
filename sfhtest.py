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
np.errstate(invalid='ignore')


if __name__ == "__main__":

    sfhtest = 0  # sfh test

    folder = 'out_sfhtest/'
    pkl = 'pkl_sfhtest/'
    # pars = ['eelg_test' + str(n) + '_params.py' for n in range(100)]

    gals = {}
    git = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
    oute = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + folder
    pkls = git + pkl
    for file in os.listdir(oute):
        c = 0
        key = ''
        for char in file:
            if char == '_':
                c += 1
            elif c == 3:
                key += char
        if file.endswith(".h5") and float(key) < 100:
            gals[key] = file

    fig = plt.figure()
    ax1 = plt.subplot2grid(shape=(2, 6), loc=(0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)
    ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)
    ax4 = plt.subplot2grid((2, 6), (1, 1), colspan=2)
    ax5 = plt.subplot2grid((2, 6), (1, 3), colspan=2)

    # get_e = np.zeros(shape=(4, len(gals)))  # 6 rows (mass, dust, logzsol, gaslogz, ssfr1, ssfr2)
    onedraw = np.zeros(shape=(len(gals), 6, 10**3))  # *10**3))  # [mass, dust, metal, gasmet, ssfr1, ssfr2]
    # NOTE: INPUT SAME gasmet each time
    # sfh_bins = np.zeros(shape=(len(gals), 2, 10**4))  # 10**4 = num arg in ssfr_fig3.randraw()
    for i in range(100):  # for each galaxy (0 through 99)
        print(i, gals[str(i)])  # i = key
        if os.path.exists(oute + gals[str(i)]):
            onedraw[i, :, :] = gmd.printer(oute + gals[str(i)], percs=False, sfhtest=True, draw1=True)

        '''
        glxy = '11063_cosmos_sfhtest_' + str(i)
        # fgal = fnames[i]
        file = pkls + glxy + '_extra_out.pkl'
        # ffile = pklsf + fgal + '_extra_out.pkl'
        randomdraws = sfh.randraw(file, onedraw[0, i, np.random.randint(10**3)])[0]  # args=(file, mass)  # size 22, num
        # note: above random randint chooses a random value of mass for galaxy

        recent = []
        second = []
        for j in range(len(randomdraws)):  # 22?
            if j == 0 or j == 1 or j == 2:
                for k in range(len(randomdraws[0])):  # num = 10**4?
                    # sfh_bins[i, 0, k] += randomdraws[j][k] / 3
                    recent.append(randomdraws[j][k])
            elif j == 3 or j == 4 or j == 5 or j == 6:
                for k in range(len(randomdraws[0])):  # num = 10**4?
                    # sfh_bins[i, 1, k] += randomdraws[j][k] / 4
                    second.append(randomdraws[j][k])
        # WAIT!!!!!!!!! COMPARING SSFR TO INPUT MFRAC ISN'T RIGHT
        '''
        pars = git + 'bettertest/sfhtest_' + str(i) + '_params.py'
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
        print(onedraw[i, 0, :])

        # ssfr050 = np.percentile(recent, [16., 50., 84.])
        # ssfr50100 = np.percentile(second, [16., 50., 84.])
        mass = np.percentile(onedraw[i, 0, :], [16., 50., 84.])
        print(mass, 'mass')
        dust = np.percentile(onedraw[i, 1, :], [16., 50., 84.])
        met = np.percentile(onedraw[i, 2, :], [16., 50., 84.])
        gasmet = np.percentile(onedraw[i, 3, :], [16., 50., 84.])
        sfr1 = np.percentile(onedraw[i, 4, :], [16., 50., 84.])
        sfr2 = np.percentile(onedraw[i, 5, :], [16., 50., 84.])

        ax1.errorbar(float(inmass), mass[1], yerr=np.array([[mass[1] - mass[0], mass[2] - mass[1]]]).T, fmt='ko')
        ax2.errorbar(float(indust), dust[1], yerr=np.array([[dust[1] - dust[0], dust[2] - dust[1]]]).T, fmt='ko')
        ax3.errorbar(10**float(inmet), 10**met[1], yerr=np.array([[10**met[1]-10**met[0], 10**met[2]-10**met[1]]]).T,
                     fmt='ko')
        ax4.errorbar(float(insfr1), sfr1[1], yerr=np.array([[sfr1[1] - sfr1[0], sfr1[2] - sfr1[1]]]).T, fmt='ko')
        ax5.errorbar(float(insfr2), sfr2[1], yerr=np.array([[sfr2[1] - sfr2[0], sfr2[2] - sfr2[1]]]).T, fmt='ko')

        # print('ssfr 0-50', ssfr050[1], '+', ssfr050[2] - ssfr050[1], '-', ssfr050[1]-ssfr050[0])
        # print('ssfr 50-100', ssfr50100[1], '+', ssfr50100[2] - ssfr50100[1], '-', ssfr50100[1]-ssfr50100[0])
        print('mass', mass[1], '+', mass[2]-mass[1], '-', mass[1]-mass[0])
        print('dust', dust[1] / 1.086, '+', (dust[2]-dust[1]) / 1.086, '-', (dust[1]-dust[0]) / 1.086)
        print('met', 10**met[1], '+', 10**met[2]-10**met[1], '-', 10**met[1]-10**met[0])
        print('gasmet', 10**gasmet[1], '+', 10**gasmet[2]-10**gasmet[1], '-', 10**gasmet[1]-10**gasmet[0])

    ax1.set_xlim(8., 10.5)
    ax1.set_ylim(8., 10.5)
    ax1.plot([8., 11.], [8., 11.], color='k', linestyle='--')
    ax2.set_xlim(0., 2.)
    ax2.set_ylim(0., 2.)
    ax2.plot([0., 2.], [0., 2.], color='k', linestyle='--')
    ax3.set_xlim(0., .2)
    ax3.set_ylim(0., .2)
    ax3.plot([0., .5], [0., .5], color='k', linestyle='--')
    ax3.set_xlim(0., .8)
    ax3.set_ylim(0., .8)
    ax4.plot([0., 1.], [0., 1.], color='k', linestyle='--')
    ax5.set_xlim(0., .5)
    ax5.set_ylim(0., .5)
    ax5.plot([0., 1.], [0., 1.], color='k', linestyle='--')
    plt.show()
