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


def density_estimation(m1, m2, xs=[8.5, 11.5], ys=[0., 700.]):
    X, Y = np.mgrid[xs[0]:xs[1]:100j, ys[0]:ys[1]:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def plotter(x, y, color, errs=True, fs_text=30, fs=20, norm=7., xs=[8.5, 11.5], ys=[0., 700.], scat=True, ax=None,
            levels=None, alpha=0.5):
    # font={'fontname': 'Times'},
    if errs:
        '''
        errx = np.percentile(x, [16., 50., 84.])
        erry = np.percentile(y, [16., 50., 84.])
        plt.errorbar(x=[errx[1]], y=[erry[1]], xerr=[[errx[1] - errx[0]], [errx[2] - errx[1]]],
                     yerr=[[erry[1] - erry[0]], [erry[2] - erry[1]]], color=color, fmt='*')
        '''
        X, Y, Z = density_estimation(x, y, xs=xs, ys=ys)
        # levels = np.arange(0.5, np.amax(Z), 0.02) + 0.02
        if levels is None:
            levels = np.arange(np.amax(Z) / norm, np.amax(Z), np.amax(Z) / norm) + (np.amax(Z) / norm)
        else:
            levels = [lv * np.amax(Z) for lv in levels]
        if color == 'b':
            cmap = 'Blues'
        else:
            cmap = 'Purples'
        if ax is not None:
            ax.contourf(X, Y, Z, levels=levels, cmap=cmap, alpha=alpha)  # [0.1, 0.2, 0.5, 1., 25.]
        else:
            plt.contourf(X, Y, Z, levels=levels, cmap=cmap, alpha=alpha)  # [0.1, 0.2, 0.5, 1., 25.]
            # plt.contour(X, Y, Z, levels=levels, cmap=cmap, lw=3)#, alpha=0.5)  # [0.1, 0.2, 0.5, 1., 25.]


if __name__ == "__main__":

    folder = 'out_sfhtest/'
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
    fig = plt.figure()
    '''
    ax1 = plt.subplot2grid(shape=(2, 6), loc=(0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)
    ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)
    ax4 = plt.subplot2grid((2, 6), (1, 1), colspan=2)
    ax5 = plt.subplot2grid((2, 6), (1, 3), colspan=2)
    '''
    ax1 = plt.subplot(2, 3, 1)
    ax2 = plt.subplot(2, 3, 2)
    ax3 = plt.subplot(2, 3, 3)
    ax4 = plt.subplot(2, 3, 4)
    ax5 = plt.subplot(2, 3, 5)
    ax6 = plt.subplot(2, 3, 6)

    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    fs_text = 20  # 30
    fs = 10  # 20
    fs_ticks = 15  # 25

    l1, l2, l3, l4, l5, l6 = [[], []], [[], []], [[], []], [[], []], [[], []], [[], []]

    # get_e = np.zeros(shape=(4, len(gals)))  # 6 rows (mass, dust, logzsol, gaslogz, ssfr1, ssfr2)
    onedraw = np.zeros(shape=(200, 6, 10**3))  # *10**3))  # [mass, dust, metal, gasmet, ssfr1, ssfr2]
    # (shape=(len(gals), 6, 10**3))
    # NOTE: INPUT SAME gasmet each time
    # sfh_bins = np.zeros(shape=(len(gals), 2, 10**4))  # 10**4 = num arg in ssfr_fig3.randraw()
    # for ind in range(len(allkeys)):
    for i in range(100):  # for each galaxy (0 through 99)
        # i = allkeys[ind]
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
        pars = git + 'ctest/sfhtest_' + str(i) + '_params.py'
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

        '''
        if 2.5 <= zred < 3.:
            color = 'b'
        elif 3. <= zred < 3.5:
            color = 'purple'
        elif 3.5 < zred < 4.:
            color = 'r'
        # '''
        color = 'k'

        ax1.errorbar(float(inmass), mass[1], yerr=np.array([[mass[1] - mass[0], mass[2] - mass[1]]]).T, fmt='o', color=color)
        ax2.errorbar(float(indust), dust[1], yerr=np.array([[dust[1] - dust[0], dust[2] - dust[1]]]).T, fmt='o', color=color)
        ax3.errorbar(10**float(inmet), 10**met[1], yerr=np.array([[10**met[1]-10**met[0], 10**met[2]-10**met[1]]]).T,
                     fmt='o', color=color)
        ax4.errorbar(float(insfr1), sfr1[1], yerr=np.array([[sfr1[1] - sfr1[0], sfr1[2] - sfr1[1]]]).T, fmt='o', color=color)
        ax5.errorbar(float(insfr2), sfr2[1], yerr=np.array([[sfr2[1] - sfr2[0], sfr2[2] - sfr2[1]]]).T, fmt='o', color=color)
        ax6.errorbar(float(insfr1) / float(insfr2), sfr1[1] / sfr2[1], fmt='o', color=color)
        # '''
        l1[0].append(float(inmass))
        l1[1].append(mass)
        l2[0].append(float(indust))
        l2[1].append(dust)
        l3[0].append(10**float(inmet))
        l3[1].append(10**met)
        l4[0].append(float(insfr1))
        l4[1].append(sfr1)
        l5[0].append(float(insfr2))
        l5[1].append(sfr2)
        l6[0].append(float(insfr1) / float(insfr2))
        l6[1].append(sfr1 / sfr2)
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
    # SET UP AXES:
    axes = [ax1, ax2, ax3, ax4, ax5]
    for ax in axes:
        ax.tick_params('x', length=3, width=1, which='both', labelsize=fs, pad=10)
        ax.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)

    # MASS AXES
    ax1.set_xlim(8., 11.)
    ax1.set_ylim(8., 11.)
    ax1.plot([8., 11.], [8., 11.], color='k', linestyle='--')
    ax1.set_xlabel(r'Input stellar mass [$\log_{10}($M/M$_\odot$)]', fontsize=fs_text)
    ax1.set_ylabel(r'Recovered stellar mass [$\log_{10}($M/M$_\odot$)]', fontsize=fs_text)
    ax1.set_xticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0])  # technically works
    ax1.set_xticklabels([r'$8.0$', r'$8.5$', r'$9.0$', r'$9.5$', r'$10.0$', r'$10.5$', r'$11.0$'], size=fs_ticks)
    ax1.set_yticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0])  # technically works
    ax1.set_yticklabels([r'$8.0$', r'$8.5$', r'$9.0$', r'$9.5$', r'$10.0$', r'$10.5$', r'$11.0$'], size=fs_ticks)

    # DUST AXES
    ax2.set_xlim(0., 2.)
    ax2.set_ylim(0., 2.)
    ax2.plot([0., 2.], [0., 2.], color='k', linestyle='--')
    ax2.set_xlabel(r'Input A$_V$', fontsize=fs_text)
    ax2.set_ylabel(r'Recovered A$_V$', fontsize=fs_text)
    ax2.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])  # technically works
    ax2.set_xticklabels([r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$'], size=fs_ticks)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])  # technically works
    ax2.set_yticklabels([r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$'], size=fs_ticks)

    # METALLICITY AXES
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlim(10**-2., 1.)
    ax3.set_ylim(10**-2., 1.)
    ax3.plot([10**-2., 1.], [10**-2., 1.], color='k', linestyle='--')
    ax3.set_xlabel(r'Input stellar mettalicity [$\log_{10}($Z/Z$_\odot$)]', fontsize=fs_text)
    ax3.set_ylabel(r'Recovered stellar mettalicity [$\log_{10}($Z/Z$_\odot$)]', fontsize=fs_text)
    ax3.set_xticks([10 ** -2, 10 ** -1, 10 ** 0])  # technically works
    ax3.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$'], size=fs_ticks)
    ax3.set_yticks([10 ** -2, 10 ** -1, 10 ** 0])  # technically works
    ax3.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$'], size=fs_ticks)

    # SFR BIN1 AXES
    ax4.set_xlim(0., 1.)
    ax4.set_ylim(0., 1.)
    ax4.plot([0., 1.], [0., 1.], color='k', linestyle='--')
    ax4.set_xlabel(r'Input fraction of M$^*$ formed [$0 - 50$ Myr]', fontsize=fs_text)
    ax4.set_ylabel(r'Recovered fraction of M$^*$ formed [$0 - 50$ Myr]', fontsize=fs_text)
    ax4.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax4.set_xticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    ax4.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax4.set_yticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)

    # SFR BIN2 AXES
    ax5.set_xlim(0., 1.)
    ax5.set_ylim(0., 1.)
    ax5.plot([0., 1.], [0., 1.], color='k', linestyle='--')
    ax5.set_xlabel(r'Input fraction of M$^*$ formed [$50 - 100$ Myr]', fontsize=fs_text)
    ax5.set_ylabel(r'Recovered fraction of M$^*$ formed [$50 - 100$ Myr]', fontsize=fs_text)
    ax5.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax5.set_xticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    ax5.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
    ax5.set_yticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)

    ax6.set_xscale("log")
    ax6.set_yscale("log")
    ax6.set_xlim(10**-2., 10**2.)
    ax6.set_ylim(10**-2., 10**2.)
    ax6.plot([10**-2., 10**2.], [10**-2., 10**2.], color='k', linestyle='--')
    ax6.set_xlabel(r'Input ratio of M$^*$ formed [$f_1 / f_2$]', fontsize=fs_text)
    ax6.set_ylabel(r'Recovered ratio of M$^*$ formed [$f_1 / f_2$]', fontsize=fs_text)
    ax6.set_xticks([10**-2., 10**-1., 10**0., 10**1., 10**2.])  # technically works
    ax6.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$'], size=fs_ticks)
    ax6.set_yticks([10**-2., 10**-1., 10**0., 10**1., 10**2.])  # technically works
    ax6.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$'], size=fs_ticks)
#    ax6.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])  # technically works
#    ax6.set_yticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$'], size=fs_ticks)
    plt.show()
