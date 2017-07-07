# -*- coding: utf-8 -*-
"""
UVJ plotter
Fluxes from rest frame catalog
Use flag from main catalog
Photometric redshifts from zout catalog
"""

import numpy as np
import matplotlib.pyplot as plt


def uvj_plot(objname, field, title=True, labels=True, lims=False, size=20):
    # UVJ plotter
    if field == 'cdfs':
        rest = '/home/jonathan/cdfs/cdfs.v1.6.9.rest.v0.9.cat'
        cat = '/home/jonathan/cdfs/cdfs.v1.6.11.cat'
        zname = '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout'
    elif field == 'cosmos':
        rest = '/home/jonathan/cosmos/cosmos.v1.3.6.rest.v0.9.cat'
        cat = '/home/jonathan/cosmos/cosmos.v1.3.8.cat'  # main catalog
        zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'  # redshift catalog
    elif field == 'uds':
        rest = '/home/jonathan/uds/uds.v1.5.8.rest.v0.9.cat'
        cat = '/home/jonathan/uds/uds.v1.5.10.cat'
        zname = '/home/jonathan/uds/uds.v1.5.8.zout'

    s = np.loadtxt(rest)  # rest frame catalog
    main = np.loadtxt(cat)  # main cat (was cdfsmaincat [now cdfsmaincat-bug])
    redshift = np.loadtxt(zname)

    objuse = []
    zcut = []
    UV = []
    VJ = []
    Ksnr = []
    table = []
    x = []
    y = []
    for i in range(len(main)):
        # USE FLAG AND Z_PHOT for each obj) (note: len(main[0] = 156; elements 153:155 = use, snr, use_nosnr, z_spec)
        objuse.append(main[i][-4])  # index good for all main cats
        if main[i][-1] >= 0:  # if z_spec exists for this galaxy
            zcut.append(main[i][-1])  # index good for all main cats
        else:
            zcut.append(redshift[i][17])  # using zout catalog, z_peak = z_phot; index good for all zout cats

        # CREATE UVJ AXES
        UV.append(-2.5 * np.log10(s[i][11] / s[i][15]))  # indices good for all rest frame cats
        VJ.append(-2.5 * np.log10(s[i][15] / s[i][17]))  # indices good for all rest frame cats
        if i == int(objname):
            special_x = -2.5 * np.log10(s[i][11] / s[i][15])
            special_y = -2.5 * np.log10(s[i][15] / s[i][17])
            print(i, objname)

        # SNR IN K BAND
        Ksnr.append(main[i][21] / main[i][22])  # f_Ksall / e_Ksall, indices good for all main cats

        # UVJ CORNER CUT, USE CUT, K-BAND SNR CUT, AND REDSHIFT CUT (zcut[i] < 2.369 if include F160W; else < 2.560)
        if objuse[i] == 1 and Ksnr[i] >= 10:  # >=20
            y.append(UV[i])
            x.append(VJ[i])
            table.append(main[i][0])  # ID in catalog = main[i][0]

    # PLOT SCATTER OF REMAINING POINTS
    plt.scatter(x, y, color='0.5', alpha=0.1, marker=".")
    plt.scatter(special_x, special_y, color='b', marker="*", s=100)

    if title:
        plt.title(field + '-' + objname)
    if labels:
        plt.text(-0.4, 1.35, 'Quiescent', fontsize=16)  # label quiescent region
        plt.text(-0.4, 1.1, 'Star-forming', fontsize=16)  # label star-forming region
    # plt.text(-0.2, 1.175, 'Star-forming', fontsize=16)  # label star-forming region
    plt.plot([-0.5, 0.9], [1.3, 1.3], 'k')  # plot a line to show quiescent region ([x1,x2], [y1, y2])
    plt.plot([1.6, 1.6], [2.5, 2.0], 'k')  # plot a line to show quiescent region ([x1,x2], [y1, y2])
    plt.plot([0.9, 1.6], [1.3, 2.0], 'k')  # plot a line to show quiescent region ([x1,x2], [y1, y2])

    if lims:
        plt.xlim(-1.5, 2.5)
        plt.xticks([-1, 0, 1, 2])
        plt.ylim(-.5, 2.5)
        plt.yticks([0.0, 1., 2.])

    plt.xlabel(r'Rest frame $V - J$', fontsize=size)
    plt.ylabel(r'Rest frame $U - V$', fontsize=size)
    plt.show()
