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
from matplotlib import rc
np.errstate(invalid='ignore')


if __name__ == "__main__":
    eels = True
    sfs = True
    path = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
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
    if eels and sfs:
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
        plt.xlabel(label, fontsize=fs_text, **font)  # 20
        plt.ylabel(ylabel, fontsize=fs_text, **font)  # 20
        plt.show()

    elif eels and not sfs:
        pass
    elif sfs and not eels:
        pass

