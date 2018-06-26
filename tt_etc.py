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

    tt = 0  # Tiantian's galaxies
    sfhtest = 1  # sfh test

    if tt:
        folders = ['out_tt/', 'out_tfast/']
        pkls = ['pkl_tt/', 'pkl_tfast/']
        pars = ['eelg_fifty_params.py', 'eelg_fast_params.py']
        base = ['tt' 'tfast']
        gnames = ['5519_cdfs_tt', '5593_cdfs_tt', '5475_cdfs_tt']
        fnames = ['5519_cdfs_tfast', '5593_cdfs_tfast', '5475_cdfs_tfast']
    elif sfhtest:
        folders = ['out_simsfh/', 'out_simsfh/']
        pkls = ['pkl_simsfh/', 'pkl_simsfh/']
        pars = ['eelg_simsfh_params.py', 'eelg_simsfh_params.py']
        base = ['9_03_17_15_05', '9_03_17_15_05']
        fnames = ['11063_sim_9_03_17_15_05']
        gnames = ['11063_sim_9']

    gals = []
    oute = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + folders[0]
    pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + pkls[0]
    for file in os.listdir(oute):
        if file.endswith(".h5"):
            gals.append(file)

    fasts = []
    outl = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + folders[1]
    pklsf = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + pkls[1]
    for file in os.listdir(outl):
        if file.endswith(".h5"):
            fasts.append(file)

    # get_e = np.zeros(shape=(4, len(gals)))  # 4 rows (dust, mass, gaslogz, logzsol), each row as long as eelgs
    onedraw = np.zeros(shape=(4, len(gals), 10**3))  # *10**3))
    onedrawf = np.zeros(shape=(4, len(fasts), 10**3))  # *10**3))
    for i in range(len(gals)):
        if os.path.exists(oute + gals[i]):
            # print('reg', printer(oute + gals[i]))
            # print('fast', printer(outl + fasts[i], percs=False, fast=True))

            onedraw[:, i, :] = gmd.printer(oute + gals[i], percs=False, draw1=True)
            onedrawf[:, i, :] = gmd.printer(outl + fasts[i], percs=False, fast=True, draw1=True)
            # for k in range(10**3):
        glxy = gnames[i]
        # fgal = fnames[i]
        file = pkls + glxy + '_extra_out.pkl'
        # ffile = pklsf + fgal + '_extra_out.pkl'
        print(i)
        randomdraws = sfh.randraw(file, onedraw[0, i, np.random.randint(10**3)])[0]  # size 22, num
        # note: above random randint chooses a random value of mass for galaxy
        recent = []
        second = []
        for j in range(len(randomdraws)):  # 22?
            if j == 0 or j == 1 or j == 2:
                for k in range(len(randomdraws[0])):  # num = 10**4?
                    recent.append(randomdraws[j][k])
            elif j == 3 or j == 4 or j == 5 or j == 6:
                for k in range(len(randomdraws[0])):  # num = 10**4?
                    second.append(randomdraws[j][k])
        print(gals[i])
        print('ssfr 0-50', np.percentile(recent, [16., 50., 84.]))
        print('ssfr 50-100', np.percentile(second, [16., 50., 84.]))
        print('mass', np.percentile(onedraw[0, i, :], [16., 50., 84.]))
        print('dust', np.percentile(onedraw[1, i, :], [16., 50., 84.]))
        print('met', np.percentile(onedraw[2, i, :], [16., 50., 84.]))
        print('gasmet', np.percentile(onedraw[3, i, :], [16., 50., 84.]))
        print(fasts[i])
        print('mass', np.percentile(onedrawf[0, i, :], [16., 50., 84.]))
        print('tage', np.percentile(onedrawf[1, i, :], [16., 50., 84.]))
        print('logtau', np.percentile(onedrawf[2, i, :], [16., 50., 84.]))
        print('dust', np.percentile(onedrawf[3, i, :], [16., 50., 84.]))

'''
5519_cdfs_tt_1529681589_mcmc.h5
('ssfr 0-50', array([  2.13811580e-09,   2.65041068e-09,   4.06467211e-09]))
('ssfr 50-100', array([  3.35884165e-10,   1.28210011e-09,   2.88002932e-09]))
('mass', array([ 8.84616123,  8.91384125,  8.98551491]))
('dust', array([ 0.09309613,  0.21738216,  0.367716  ]))
('met', array([-0.45945215, -0.12037901,  0.09563315]))
('gasmet', array([-1.57574482, -0.74148431, -0.1373779 ]))
5519_cdfs_tfast_1529684132_mcmc.h5
('mass', array([ 8.80425129,  8.85554743,  8.89594406]))
('tage', array([ 2.45612302,  3.50806594,  4.18618896]))
('logtau', array([-1.34073757,  0.10544156,  1.38699596]))
('dust', array([ 0.32267835,  0.37998955,  0.47563315]))

5593_cdfs_tt_1529686734_mcmc.h5
('ssfr 0-50', array([  6.97110150e-10,   8.26257842e-10,   9.83950533e-10]))
('ssfr 50-100', array([  5.47038084e-11,   2.46552797e-10,   7.36842120e-10]))
('mass', array([ 8.10925735,  8.21538925,  8.30441715]))
('dust', array([ 0.03640653,  0.10723628,  0.20377764]))
('met', array([-0.0227099 ,  0.09697676,  0.161579  ]))
('gasmet', array([-1.80078944, -1.40309417, -1.1121323 ]))
5593_cdfs_tfast_1529689965_mcmc.h5
('mass', array([ 8.04089443,  8.09354305,  8.14443062]))
('tage', array([ 2.99748585,  3.92401087,  4.79658571]))
('logtau', array([-1.45291759, -0.17301916,  1.31309599]))
('dust', array([ 0.033387  ,  0.10405507,  0.19853153]))

5475_cdfs_tt_1529681672_mcmc.h5
('ssfr 0-50', array([  7.65770823e-11,   1.17383962e-10,   1.70387840e-10]))
('ssfr 50-100', array([  1.37224271e-11,   6.51601882e-11,   1.95537203e-10]))
('mass', array([ 9.51058006,  9.58273792,  9.64859978]))
('dust', array([ 0.03305917,  0.12124067,  0.26184136]))
('met', array([-0.23366278, -0.00958648,  0.13142096]))
('gasmet', array([-0.89538906, -0.08187881,  0.28848711]))
5475_cdfs_tfast_1529684086_mcmc.h5
('mass', array([ 9.38911087,  9.43760014,  9.48172237]))
('tage', array([ 1.67573083,  2.15706217,  2.42894878]))
('logtau', array([-1.28101252,  0.04709384,  1.38593617]))
('dust', array([ 0.31577208,  0.36668095,  0.44947072]))
'''