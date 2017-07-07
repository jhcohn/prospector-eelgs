import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import print_sfh
import uvj
from matplotlib import gridspec


def all_plots(files, objname, field, loc='upper left'):
    # PLOT UVJ
    # uvj.uvj_plot(objname, field)

    # PLOT SFH, SED+
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
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])  # prepares fig to hold two axes, one atop other, size ratio 3:1

    ax1 = plt.subplot(gs[0])  # ax1 = bigger upper axis
    # plt.axvline(x=5270, color='k')  # proving no offset in two axes
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_title(field + '-' + obj)
    line11 = ax1.errorbar(wave_rest, results['obs']['maggies'], yerr=results['obs']['maggies_unc'],
                          marker='o', linestyle='', color='b')  # plot observations
    line12 = ax1.plot(wave_rest, sed, 'o', color='r')  # plot best fit model
    line13 = ax1.plot(sps_wave, spec, color='b', alpha=0.5)  # plot spectrum
    ax1.set_ylabel(r'Maggies')

    ax2 = plt.subplot(gs[1], sharex=ax1)  # ax2 = smaller lower axis
    line2 = ax2.plot(wave_rest, chi, 'o', color='k')  # plot chi
    plt.axhline(y=0, color='k')  # plot horizontal line at y=0 on chi plot
    plt.setp(ax1.get_xticklabels(), visible=False)  # hide xtick labels on upper axis
    yticks = ax2.yaxis.get_major_ticks()  # show ytick labels on lower axis
    yticks[-1].label1.set_visible(False)  # hide uppermost ytick label on lower axis to prevent overlap
    ax2.set_ylabel(r'$\chi$')
    ax2.set_xlabel(r'Rest frame wavelength')

    ax1.legend((line11, line12, line13, line2), loc=loc, prop={'size': 14},
               labels=[r'Model', r'Spectrum', r'Observed Photometry', r'$\chi$'])

    # plt.axvline(x=5270, color='k')  # proving no offset in two axes
    plt.subplots_adjust(hspace=.0)

    # TESTING
    # UVJ inset on SED+ plot
    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    inset = inset_axes(ax1, width="25%", height=2.5, loc=2)  # create inset axis: width (%), height (inches), location
    uvj.uvj_plot(objname, field, title=False, labels=False, lims=True, size=20)  # add uvj plot to inset axis
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

    all_plots(files, obj, field, loc='upper right')

'''
Currently running with:
python make_all_plots.py --obj=17423 --field=cosmos
'''

'''
plt.subplot(211, xscale="log", yscale="log")
plt.errorbar(wave_rest, results['obs']['maggies'], yerr=results['obs']['maggies_unc'],
             marker='o', linestyle='', color='b', label='Observed photometry')
plt.plot(wave_rest, sed, 'o', color='r')
plt.plot(sps_wave, spec, color='b', alpha=0.5)

plt.subplot(212)
plt.plot(wave_rest, chi_sq, 'o', color='b')

#                            plt.errorbar(wave_rest, results['obs']['maggies'], yerr=results['obs']['maggies_unc'],
#                                         marker='o', linestyle='', color='b', label='Observed photometry')
#                            plt.plot(wave_rest, sed, 'o', label='Model at {},{}'.format(walker, iteration), color='r')
#                            plt.legend(loc="best", fontsize=20)
#                            plt.title(str(objname) + ' SED')
#                            plt.plot(sps_wave, spec, color='b', alpha=0.5)
#                            plt.xlabel('Rest frame wavelength [angstroms]')
#                            plt.ylabel('Maggies')
plt.show()
'''
'''
plt.plot(wave_rest, chi_sq, 'o', color='b')
plt.title(str(objname) + r' $\chi^2$')
plt.xlabel('Rest frame wavelength [angstroms]')
plt.ylabel(r'$\chi^2$')
plt.show()
'''