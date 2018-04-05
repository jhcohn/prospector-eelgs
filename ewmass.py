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
import glob
import copy
from astropy import constants
import eelg_fifty_params as pfile
import get_mass_dust as gmd
import stellar_ages as sa
from matplotlib import rc
from scipy import stats
np.errstate(invalid='ignore')
# NOTE: check if OIII / OII high for SFGs because OII essentially just noisy around 0?


def list_ssfr(file):
    '''
    :param file: .txt file listing galaxy id followed by the galaxy's ssfr
    :return:
    '''
    galaxies = []
    counter = 0
    with open(file, 'r') as ssfr:
        for line in ssfr:
            if line[0] != '#':
                cols = line.split()
                # print(cols[0])
                galaxies.append(float(cols[1]) * 10 ** 9)  # ssfrs are stored in 2nd column in file; convert to Gyr^-1
                counter += 1
    return galaxies


def list_ewflux(file, line_idx, ratio=False):
    '''
    lines in files are listed in order: '[NII]6585', '[OII]3726', r'H$\alpha$6563', '[OII]3729', '[OIII]4960',
    '[OIII]5007', r'H$\beta$4861',  r'H$\delta$4102'

    :param file: .txt file listing galaxy id followed by the galaxy's emission line widths OR line fluxes
    :param line_idx: idx 0 - 7 for the desired line EW or flux, OR [idx1, idx2] of lines for desired flux ratio
    :return:
    '''
    galaxies = []
    counter = 0
    with open(file, 'r') as ewflux:
        for line in ewflux:
            if line[0] != '#':
                cols = line.split()
                # print(cols[0])
                if not ratio:
                    galaxies.append(float(cols[line_idx + 1]))  # 0th line_idx is ID
                else:
                    print(line_idx, line_idx[0])
                    print(cols)
                    galaxies.append(float(cols[line_idx[0] + 1]) / float(cols[line_idx[1] + 1]))
                counter += 1

    return galaxies


def preplot(order=None, keys='e', file=None, line_idx=5, line_idx2=[5, 3], ratio=[False, False], outfold=None,
            fs_text=30, fs=20, color='purple', font={'fontname': 'Times'}):
    '''

    :param keys:
    :param files:
    :param line_idx: single idx for desired line flux or EW, OR [idx1, idx2] for desired line ratios
    :param line_idx2: (if no ssfr) single idx for desired line flux or EW, OR [idx1, idx2] for desired line ratios
    :param ratio:
    :param fs_text:
    :param fs:
    :param font:
    :return:
    '''

    lines = ['[NII]6585', '[OII]3726', r'H$\alpha$6563', '[OII]3729', '[OIII]4960', '[OIII]5007', r'H$\beta$4861',
             r'H$\delta$4102']
    # names = '[NII]6585', '[OII]3726', 'Halpha6563', '[OII]3729', '[OIII]4960', '[OIII]5007', 'Hbeta4861', 'Hdelta4102'

    labels = [None, None]
    ew = []
    if keys != 's':
        ew = list_ewflux(file, line_idx=line_idx, ratio=ratio[0])  # line1
        # xy[1] = list_ewflux(files[1], line_idx=line_idx2, ratio=ratio[1])  # line2
        labels = [r'Stellar Mass [$\log_{10}(\textrm{M} / \textrm{M}_\odot$)]',
                  lines[line_idx] + r' Equivalent Width [\AA]']
        # labels[1] = lines[line_idx2[0]] + '/' + lines[line_idx2[1]]
    else:  # want ssfr
        ew = list_ssfr(file)  # ssfr
        labels = [r'Stellar Mass [log_{10}(M / M$_\odot$)]', 'SSFR [Gyr$^{-1}$]']

    out = outfold
    gals1 = []
    for file in os.listdir(out):
        if file.endswith(".h5"):
            gals1.append(file)
    print(gals1, order)
    gals = []
    for ln in range(len(order)):
        for id in range(len(gals1)):
            if gals1[id].startswith(str(order[ln])):
                gals.append(gals1[id])
    print(gals, order)

    mass = []
    for i in range(len(gals1)):
        mass.append(gmd.printer(out + gals1[i], percs=True)[0])  # median mass

    return ew, mass, labels


def density_estimation(m1, m2):
    X, Y = np.mgrid[8.5:11.5:100j, 0.:700.:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def plotter(x, y, color, fs_text=30, fs=20, font={'fontname': 'Times'}):

    plt.scatter(x=x, y=y, color=color)

    errx = np.percentile(x, [16., 50., 84.])
    erry = np.percentile(y, [16., 50., 84.])
    plt.errorbar(x=[errx[1]], y=[erry[1]], xerr=[[errx[1] - errx[0]], [errx[2] - errx[1]]],
                 yerr=[[erry[1] - erry[0]], [erry[2] - erry[1]]], color=color, fmt='*')

    X, Y, Z = density_estimation(x, y)
    # levels = np.arange(0.5, np.amax(Z), 0.02) + 0.02
    levels = np.arange(np.amax(Z) / 100., np.amax(Z), np.amax(Z) / 10.) + (np.amax(Z) / 10.)
    if color == 'b':
        cmap = 'Blues'
    else:
        cmap = 'Purples'
    plt.contourf(X, Y, Z, levels=levels, cmap=cmap, alpha=0.5)  # levels=levels  # [0.1, 0.2, 0.5, 1., 25.]


if __name__ == "__main__":

    path = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
    e_out = path + 'out/out_efico/'
    s_out = path + 'out/out_nfico/'
    e_fs = {'s': path + 'eelg_ssfrs.txt', 'e': path + 'ews_eelgs.txt', 'f': path + 'fluxes_eelgs.txt'}
    s_fs = {'s': path + 'sfg_ssfrs.txt', 'e': path + 'ews_sfgs.txt', 'f': path + 'fluxes_sfgs.txt'}

    # same for EELGs, SFGs
    key = 'e'
    ratio = [False, True]
    idx = 5  # [NII]6585, [OII]3726, 'Halpha6563', '[OII]3729', '[OIII]4960', '[OIII]5007', 'Hbeta4861', 'Hdelta4102'

    ee_order, ee_fd, sf_order, sf_fd = sa.get_gal_lists(base=['fico', 'fico'], objlists=True, normal=True)

    # different for EELGs, SFGs
    colors = ['purple', 'b']  # eelg, sfg
    eelg_file = e_fs[key]
    sfg_file = s_fs[key]

    eew, emass, labs = preplot(order=ee_order, keys=key, file=eelg_file, line_idx=idx, ratio=ratio, outfold=e_out)
    print('hi')
    sew, smass, labs = preplot(order=sf_order, keys=key, file=sfg_file, line_idx=idx, ratio=ratio, outfold=s_out)
    print(len(eew), len(emass), len(sew), len(smass))

    sfg2_x = []
    sfg2_y = []
    for thing in range(len(smass)):
        if sew[thing] <= 10**-16 or smass[thing] <= 10**-16:
            # print('me')
            pass
        else:
            sfg2_x.append(smass[thing])
            sfg2_y.append(sew[thing])
    smass = sfg2_x
    sew = sfg2_y
    print(len(eew), len(emass), len(sew), len(smass))

    plotter(emass, eew, color=colors[0])
    plotter(smass, sew, color=colors[1])

    # PLOT PARAMS
    fs_text = 30
    fs = 20
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    font = {'fontname': 'Times'}
    xmin, xmax = 8.5, 11.5  # 1.2 * max([max(emass), max(smass)])  # 50.  # 10 ** 3  # 700.
    ymin, ymax = 0., 1.2 * max([max(eew), max(sew)])  # 25.
    log = False
    if log:
        plt.yscale('log')
        plt.xscale('log')
        xmin, xmax = 10**-2, 10**3
        ymin, ymax = 3*10**-2, 3*10
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    plt.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    plt.xlabel(labs[0], fontsize=fs_text, **font)  # 20
    plt.ylabel(labs[1], fontsize=fs_text, **font)  # 20
    plt.show()

    '''
    line_idx = 6  # Hbeta  # 2  # Halpha  # 5  # OIII_5007
    lines = ['[NII]6585', '[OII]3726', r'H$\alpha$6563', '[OII]3729', '[OIII]4960', '[OIII]5007', r'H$\beta$4861',
             r'H$\delta$4102']
    # names = '[NII]6585', '[OII]3726', 'Halpha6563', '[OII]3729', '[OIII]4960', '[OIII]5007', 'Hbeta4861', 'Hdelta4102'

    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    font = {'fontname': 'Times'}
    fs_text = 30
    fs = 20

    eelgs = np.zeros(shape=(19, 2))  # 19 galaxies, ews & ssfrs
    sfgs = np.zeros(shape=(167, 2))
    r32_e = []  # [OIII] / [OII]
    r32_s = []
    rOb_e = []  # [OIII]5007 / Hbeta
    rOb_s = []
    counter = 0
    counter2 = 0
    with open(path + 'ews_eelgs.txt', 'r') as eelg_ews:
        for line in eelg_ews:
            if line[0] != '#':
                cols = line.split()
                # print(cols[0])
                eelgs[counter, 0] = float(cols[line_idx + 1])  # 0th line_idx is ID
                counter += 1
                r32_e.append((float(cols[5+1]) + float(cols[4+1])) / (float(cols[3+1]) + float(cols[1+1])))
                rOb_e.append(float(cols[5 + 1]) / float(cols[6 + 1]))
    with open(path + 'eelg_ssfrs.txt', 'r') as eelg_ssfr:
        counter = 0
        for line in eelg_ssfr:
            if line[0] != '#':
                cols = line.split()
                # print(cols[0])
                eelgs[counter, 1] = float(cols[1])  # ssfrs are stored in second column
                counter += 1
    print(eelgs)

    with open(path + 'ews_sfgs.txt', 'r') as sfg_ews:
        for line in sfg_ews:
            if line[0] != '#':
                cols = line.split()
                # print(cols[0])
                sfgs[counter2, 0] = float(cols[line_idx + 1])  # 0th line_idx is ID
                counter2 += 1
                r32_s.append((float(cols[5+1]) + float(cols[4+1])) / (float(cols[3+1]) + float(cols[1+1])))
                rOb_s.append(float(cols[5 + 1]) / float(cols[6 + 1]))
    with open(path + 'sfg_ssfrs.txt', 'r') as sfg_ssfr:
        counter2 = 0
        for line in sfg_ssfr:
            if line[0] != '#':
                cols = line.split()
                sfgs[counter2, 1] = float(cols[1])  # ssfrs are stored in second column
                counter2 += 1
    print(sfgs)

    sfgs2 = np.zeros(shape=(127, 2))
    r32_s2 = []
    rOb_s2 = []
    for i in range(len(sfgs)):
        if sfgs[i, 1] <= 10**-16:
            # print('me')
            pass
        else:
            sfgs2[i, :] = sfgs[i, :]
            r32_s2.append(r32_s[i])
            rOb_s2.append(rOb_s[i])
    sfgs = sfgs2
    r32_s = np.asarray(r32_s2)
    rOb_s = rOb_s2
    print(sfgs)

    r32 = False
    Ob = True
    fig = plt.figure()
    x_e = eelgs[:, 0]
    y_e = eelgs[:, 1] * 10 ** 9
    x_s = sfgs[:, 0]
    y_s = sfgs[:, 1] * 10 ** 9
    cut = False
    if cut:  # gets messy below EWs of 100 Angstroms, so fit line only to data > 100 Angstroms
        x = []
        y = []
        for i in range(len(x_e)):
            if x_e[i] > 100.:
                x.append(x_e[i])
                y.append(y_e[i])
        x2 = []
        y2 = []
        for i in range(len(x_s)):
            if x_s[i] > 100.:
                x2.append(x_s[i])
                y2.append(y_s[i])
        plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), color='purple', linestyle=':')

    if r32:
        print(len(r32_e), len(y_e), len(r32_s), len(y_s))
        plt.scatter(x=r32_e, y=y_e, color='purple')
        plt.scatter(x=r32_s, y=y_s, color='b')

        e_r = np.percentile(r32_e, [16., 50., 84.])
        e_ssfr = np.percentile(y_e, [16., 50., 84.])
        plt.errorbar(x=[e_r[1]], y=[e_ssfr[1]], xerr=[[e_r[1] - e_r[0]], [e_r[2] - e_r[1]]],
                     yerr=[[e_ssfr[1] - e_ssfr[0]], [e_ssfr[2] - e_ssfr[1]]], color='purple', fmt='*')

        s_r = np.percentile(r32_s, [16., 50., 84.])
        s_ssfr = np.percentile(y_s, [16., 50., 84.])
        plt.errorbar(x=[s_r[1]], y=[s_ssfr[1]], xerr=[[s_r[1] - s_r[0]], [s_r[2] - s_r[1]]],
                     yerr=[[s_ssfr[1] - s_ssfr[0]], [s_ssfr[2] - s_ssfr[1]]], color='b', fmt='*')
        xmin, xmax = 0., 50.  # 10 ** 3  # 700.
        ymin, ymax = 0., 25.
        ylabel = 'SSFR [Gyr$^{-1}$]'
        label = r'[OIII] / [OII]'
        # label = lines[line_idx] + ' Equivalent Width [\AA]'
    elif Ob:
        print(len(r32_e), len(y_e), len(r32_s), len(y_s))
        plt.scatter(x=x_e, y=rOb_e, color='purple')
        plt.scatter(x=x_s, y=rOb_s, color='b')

        e_r = np.percentile(x_e, [16., 50., 84.])
        e_ssfr = np.percentile(rOb_e, [16., 50., 84.])
        plt.errorbar(x=[e_r[1]], y=[e_ssfr[1]], xerr=[[e_r[1] - e_r[0]], [e_r[2] - e_r[1]]],
                     yerr=[[e_ssfr[1] - e_ssfr[0]], [e_ssfr[2] - e_ssfr[1]]], color='purple', fmt='*')

        s_r = np.percentile(x_s, [16., 50., 84.])
        s_ssfr = np.percentile(rOb_s, [16., 50., 84.])
        plt.errorbar(x=[s_r[1]], y=[s_ssfr[1]], xerr=[[s_r[1] - s_r[0]], [s_r[2] - s_r[1]]],
                     yerr=[[s_ssfr[1] - s_ssfr[0]], [s_ssfr[2] - s_ssfr[1]]], color='b', fmt='*')
        xmin, xmax = 0., 150.  # 10 ** 3  # 700.
        ymin, ymax = 0., 6.
        ylabel = r'[OIII] / H$\beta$'
        label = lines[line_idx] + ' Equivalent Width [\AA]'
    else:
        plt.plot(np.unique(x_e), np.poly1d(np.polyfit(x_e, y_e, 1))(np.unique(x_e)), color='purple', linestyle='--')
        # plt.plot(np.unique(x2), np.poly1d(np.polyfit(x2, y2, 1))(np.unique(x2)), color='b', linestyle=':')
        plt.plot(np.unique(x_s), np.poly1d(np.polyfit(x_s, y_s, 1))(np.unique(x_s)), color='b', linestyle='--')

        plt.scatter(x=x_e, y=y_e, color='purple')  # x=EW, y=ssfrs
        plt.scatter(x=x_s, y=y_s, color='b')
        e_ew = np.percentile(x_e, [16., 50., 84.])
        e_ssfr = np.percentile(y_e, [16., 50., 84.])
        plt.errorbar(x=[e_ew[1]], y=[e_ssfr[1]], xerr=[[e_ew[1] - e_ew[0]], [e_ew[2] - e_ew[1]]],
                     yerr=[[e_ssfr[1] - e_ssfr[0]], [e_ssfr[2] - e_ssfr[1]]], color='purple', fmt='*')
        s_ew = np.percentile(x_s, [16., 50., 84.])
        # print(s_ew)
        s_ssfr = np.percentile(y_s, [16., 50., 84.])
        # print(s_ssfr)
        plt.errorbar(x=[s_ew[1]], y=[s_ssfr[1]], xerr=[[s_ew[1] - s_ew[0]], [s_ew[2] - s_ew[1]]],
                     yerr=[[s_ssfr[1] - s_ssfr[0]], [s_ssfr[2] - s_ssfr[1]]], color='b', fmt='*')
        xmin, xmax = 0., 200.  # 200 [Hbeta]  # 10**3 [Ha]  # 700. [OIII]
        ymin, ymax = 0., 20.
        ylabel = 'SSFR [Gyr$^{-1}$]'
        label = lines[line_idx] + ' Equivalent Width [\AA]'
    '''
