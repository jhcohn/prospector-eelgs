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
