import eelg_multirun_params
import numpy as np
from sedpy import observate
import matplotlib.pyplot as plt


def plot_filts(ax, field, zred, scale=1, rest=True):
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
                    # y.append(float(line.split()[2]) * scale)  # scale the height of curves, default normalized to 1
                    y.append(float(line.split()[2]))  # scale the height of the curves, default normalized to 1
            y = [0.2 * y[i] / max(y) for i in range(len(y))]
            ax.fill(x, y, alpha=0.2)
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
    for thing in ids:
        new.append(filterset[thing])
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

    # vertical lines to mark OIII+H-beta region
    plt.axvspan(4800, 5050, color='k', alpha=0.3)
    plt.xlabel(r'Rest frame wavelength')
    plt.ylabel(r'Maggies')
    # plt.show()


def fig2(ax, field, zred, scale=1, rest=True):
    filterset = []
    base = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/sedpy-0.1.0-py2.7.egg/sedpy/data/filters/'

    if field == 'cdfs':
        for name in eelg_multirun_params.cdfs_filts:
            filterset.append(base + name + '.par')
            ids = [0, 4, 1, 5, 8, 9, 10, 6, 7, 11, -4, -3, -2, -1, 14]  # BVIZ, Hs Hl J1 J2 J3 Ks, IRACx4, Nb209
            colors = ['y', 'b', 'g', 'r', 'k', 'b', 'g', 'm', 'y', 'r', 'c', 'm', 'y', 'k', 'b']
    elif field == 'cosmos':
        for name in eelg_multirun_params.cos_filts:
            filterset.append(base + name + '.par')
            ids = [11, 1, 10, 2, 18, 19, 20, 17, 16, 21, -4, -3, -2, -1, 23]  # UGRI, Hs Hl J1 J2 J3 Ks, IRACx4, Nb209
            colors = ['k', '#e67e22', 'm', 'c', 'k', 'b', 'g', 'm', 'y', 'r', 'c', 'm', 'y', 'k', 'b']
            # orange
    # ids = [1, 2, 10, 11, 14, 16, 17, 18, 19, 20, 21, 22, 23, -4, -3, -2, -1]  # GIRUZ, Hl Hs J1 J2 J3 Ks Nbx2, IRACx4
    elif field == 'uds':
        for name in eelg_multirun_params.uds_filts:
            filterset.append(base + name + '.par')
            ids = [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, -4, -3, -2, -1]  # BRiz, Hs Hl J1 J2 J3 Ks, IRACx4
            # colors = ['c', 'm', 'k', 'b', 'k', 'b', 'g', 'm', 'y', 'r', 'c', 'm', 'y', 'k', 'b']
            colors = ['g', 'r', 'k', 'y', 'k', 'b', 'g', 'm', 'y', 'r', 'c', 'm', 'y', 'k', 'b']
    new = []
    for thing in ids:
        new.append(filterset[thing])
    new = np.array(new)
    # ['b', 'r', 'g', 'purple', 'k', 'purple', 'g', 'pink', 'yellow', 'r', 'purple', 'b', 'pink', 'yellow', 'k']

    i = 0
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
            ax.fill(x, y, colors[i], alpha=0.2)  # .25
            i += 1
