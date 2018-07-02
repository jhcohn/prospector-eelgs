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
    ttring = 1
    sfhtest = 0  # sfh test

    if tt:
        folders = ['out_tt/', 'out_tfastnoem/']
        pkls = ['pkl_tt/', 'pkl_tfn/']
        pars = ['eelg_fifty_params.py', 'eelg_fastnoem_params.py']
        base = ['tt' 'tfn']
        gnames = ['5519_cosmos_tt', '5593_cosmos_tt', '5475_cosmos_tt']
        fnames = ['5519_cosmos_tfn', '5593_cosmos_tfn', '5475_cosmos_tfn']

    elif ttring:
        folders = ['out_ttring2/', 'out_tfastnoem/']
        pkls = ['pkl_ttring2/', 'pkl_tfn/']
        pars = ['eelg_ttring_params.py', 'eelg_fastnoem_params.py']
        base = ['ttring2' 'tfn']
        gnames = ['5519001_cosmos_ttring2', '5519002_cosmos_ttring2', '5519003_cosmos_ttring2', '5519_cosmos_ttring2']
        fnames = ['5519_cosmos_tfn', '5593_cosmos_tfn', '5475_cosmos_tfn', '5519_cosmos_tfn']

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

        ssfr050 = np.percentile(recent, [16., 50., 84.])
        ssfr50100 = np.percentile(second, [16., 50., 84.])
        mass = np.percentile(onedraw[0, i, :], [16., 50., 84.])
        dust = np.percentile(onedraw[1, i, :], [16., 50., 84.])
        met = np.percentile(onedraw[2, i, :], [16., 50., 84.])
        gasmet = np.percentile(onedraw[3, i, :], [16., 50., 84.])

        print('ssfr 0-50', ssfr050[1], '+', ssfr050[2] - ssfr050[1], '-', ssfr050[1]-ssfr050[0])
        print('ssfr 50-100', ssfr50100[1], '+', ssfr50100[2] - ssfr50100[1], '-', ssfr50100[1]-ssfr50100[0])
        print('mass', mass[1], '+', mass[2]-mass[1], '-', mass[1]-mass[0])
        print('dust', dust[1] / 1.086, '+', (dust[2]-dust[1]) / 1.086, '-', (dust[1]-dust[0]) / 1.086)
        print('met', 10**met[1], '+', 10**met[2]-10**met[1], '-', 10**met[1]-10**met[0])
        print('gasmet', 10**gasmet[1], '+', 10**gasmet[2]-10**gasmet[1], '-', 10**gasmet[1]-10**gasmet[0])
        '''
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