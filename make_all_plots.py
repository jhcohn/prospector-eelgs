import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import print_sfh
from matplotlib import gridspec


def all_plots(files):
    # files should be a list of output files in order: sfh, res, sed, restwave, spec, spswave, chisq, justchi
    print_sfh.plotter(files[0])

    with open(files[1], 'rb') as res:
        results = pickle.load(res)

    with open(files[2], 'rb') as sed:
        sed = pickle.load(sed)

    with open(files[3], 'rb') as restwave:
        wave_rest = pickle.load(restwave)

    with open(files[4], 'rb') as spec:
        spec = pickle.load(spec)

    with open(files[5], 'rb') as spswave:
        sps_wave = pickle.load(spswave)

    '''
    with open(files[6], 'rb') as chisq:
        chi_sq = pickle.load(chisq)
    '''
    with open(files[7], 'rb') as justchi:
        chi = pickle.load(justchi)

    # plt.subplot(111, xscale="log", yscale="log")
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])

    ax1 = plt.subplot(gs[0])
    # plt.axvline(x=5270, color='k')
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_title(field + '-' + obj)
    line11 = ax1.errorbar(wave_rest, results['obs']['maggies'], yerr=results['obs']['maggies_unc'],
                          marker='o', linestyle='', color='b')  # , label='Observed photometry')
    line12 = ax1.plot(wave_rest, sed, 'o', color='r')
    line13 = ax1.plot(sps_wave, spec, color='b', alpha=0.5)
    ax1.set_ylabel('Maggies')

    ax2 = plt.subplot(gs[1], sharex=ax1)
    line2 = ax2.plot(wave_rest, chi, 'o', color='k')  # chi_sq
    zero = [0] * len(sps_wave)
    ax2.plot(sps_wave, zero, color='k')
    plt.setp(ax1.get_xticklabels(), visible=False)
    yticks = ax2.yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    # ax2.set_ylabel(r'$\chi^2$')
    ax2.set_ylabel(r'$\chi$')
    ax2.set_xlabel('Rest frame wavelength')

    ax1.legend((line11, line12, line13, line2), loc='upper right',
               labels=['Model', 'Spectrum', 'Observed Photometry', r'$\chi$'])

    # plt.axvline(x=5270, color='k')
    plt.subplots_adjust(hspace=.0)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--obj')
    parser.add_argument('--field')

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    obj = kwargs['obj']

    field = kwargs['field']
    pre = obj + '_' + field
    base = '_out.pkl'

    sfh = pre + '_sfh' + base
    res = pre + '_res' + base
    sed = pre + '_sed' + base
    restwave = pre + '_restwave' + base
    spec = pre + '_spec' + base
    spswave = pre + '_spswave' + base
    chisq = pre + '_chisq' + base
    justchi = pre + '_justchi' + base

    files = [sfh, res, sed, restwave, spec, spswave, chisq, justchi]

    all_plots(files)

'''
Currently running with:
python make_all_plots.py --obj=17423 --field=cosmos
'''