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
np.errstate(invalid='ignore')


def print_kl(out_file, true_obj):
    # ''' #
    res, pr, mod = bread.results_from(out_file)

    # HOW CONVERGED IS THE CODE?? LET'S FIND OUT!
    parnames = np.array(res['model'].theta_labels())
    fig, kl_ax = plt.subplots(1, 1, figsize=(7, 7))
    for i in xrange(parnames.shape[0]):
        kl_ax.plot(res['kl_iteration'], np.log10(res['kl_divergence'][:, i]),
                   'o', label=parnames[i], lw=1.5, linestyle='-', alpha=0.6)

    kl_ax.set_ylabel('log(KL divergence)')
    kl_ax.set_xlabel('iteration')
    # kl_ax.set_xlim(0, nsteps*1.1)

    kl_div_lim = res['run_params'].get('convergence_kl_threshold', 0.018)
    kl_ax.axhline(np.log10(kl_div_lim), linestyle='--', color='red', lw=2, zorder=1)

    kl_ax.legend(prop={'size': 5}, ncol=2, numpoints=1, markerscale=0.7)
    plt.title(true_obj + ' kl')
    plt.show()
    # ''' #


if __name__ == "__main__":
    # don't create keyword if not passed in!

    with open('Comp_10.dat', 'r') as comp:
        e_objs = []
        e_fs = []
        for line in comp:
            if line[0] == '#':
                pass
            else:
                cols = line.split()
                if int(cols[0]) - 200000 > 0:
                    e_objs.append(str(int(cols[0]) - 200000))
                    e_fs.append('uds')
                elif int(cols[0]) - 100000 > 0:
                    e_objs.append(str(int(cols[0]) - 100000))
                    e_fs.append('cosmos')
                else:
                    e_objs.append(str(int(cols[0])))
                    e_fs.append('cdfs')
    print(e_objs)
    for i in range(len(e_objs)):
        print(e_objs[i])
        for infile in glob.glob(os.path.join('/home/jonathan/.conda/envs/snowflakes/lib/python2.7/' +
                                             'site-packages/prospector/git/out/out_evar', e_objs[i] + '*.h5')):
            print(infile)
            true_field = e_fs[i]
            true_obj = e_objs[i]

            print_kl(infile, true_obj)

'''
RUNNING WITH:
python quickgrab.py --outname=10246_cdfs_multirun_1498677216_mcmc.h5 --parfile=eelg_multirun_params.py
'''