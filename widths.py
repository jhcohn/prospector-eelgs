import eelg_multirun_params
import numpy as np
from sedpy import observate
import matplotlib.pyplot as plt


def plot_filts(field, zred, scale=1, rest=True):
    filterset = []
    base = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/sedpy-0.1.0-py2.7.egg/sedpy/data/filters/'

    if field == 'cosmos':
        for name in eelg_multirun_params.cos_filts:
            filterset.append(base + name + '.par')
    elif field == 'cdfs':
        for name in eelg_multirun_params.cdfs_filts:
            filterset.append(base + name + '.par')
    elif field == 'uds':
        for name in eelg_multirun_params.uds_filts:
            filterset.append(base + name + '.par')

    filterset = np.array(filterset)

    for curve in filterset:  # for each filter in filterset
        with open(curve, 'r') as file:
            x = []
            y = []
            for line in file:
                if line[0] == 'K':  # only care about lines starting with 'KFILTER' (this skips over comments etc.)
                    if rest:  # if want rest-frame instead of observed frame
                        x.append(float(line.split()[1]) / (1 + zred))  # convert from observed to rest frame
                    else:
                        x.append(float(line.split()[1]))
                    y.append(float(line.split()[2]) * scale)  # scale the height of the curves, default normalized to 1
            plt.fill(x, y, alpha=0.25)
    # plt.show()


def some_filts(field, zred, scale=1, rest=True):
    filterset = []
    base = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/sedpy-0.1.0-py2.7.egg/sedpy/data/filters/'

    if field == 'cosmos':
        for name in eelg_multirun_params.cos_filts:
            filterset.append(base + name + '.par')
            ids = [1, 2, 10, 11, 16, 17, 18, 19, 20, 21, 23, -4, -3, -2, -1]  # GIRU, Hl Hs J1 J2 J3 Ks Nb209, IRACx4
    # ids = [1, 2, 10, 11, 14, 16, 17, 18, 19, 20, 21, 22, 23, -4, -3, -2, -1]  # GIRUZ, Hl Hs J1 J2 J3 Ks Nbx2, IRACx4
    elif field == 'cdfs':
        for name in eelg_multirun_params.cdfs_filts:
            filterset.append(base + name + '.par')
            ids = [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 14, -4, -3, -2, -1]  # BIVZ, Hs Hl J1 J2 J3 Ks Nb209, IRACx4
    elif field == 'uds':
        for name in eelg_multirun_params.uds_filts:
            filterset.append(base + name + '.par')
            ids = [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, -4, -3, -2, -1]  # BRiz, J1 J2 J3 Hs Hl Ks, IRACx4

    new = []
    for id in ids:
        new.append(filterset[id])
    new = np.array(new)

    for curve in new:  # for each filter in filterset
        with open(curve, 'r') as file:
            x = []
            y = []
            for line in file:
                if line[0] == 'K':  # only care about lines starting with 'KFILTER' (this skips over comments etc.)
                    if rest:  # if want rest-frame instead of observed frame
                        x.append(float(line.split()[1]) / (1 + zred))  # convert from observed to rest frame
                    else:
                        x.append(float(line.split()[1]))
                    y.append(float(line.split()[2]) * scale)  # scale the height of the curves, default normalized to 1
            plt.fill(x, y, alpha=0.25)
    plt.xlabel(r'Rest frame wavelength')
    plt.ylabel(r'Maggies')
    # plt.show()

'''
# EDIT to select which specific filters I want to show
'''