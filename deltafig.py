import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
import argparse
import print_sfh
import uvj
import widths  # WIDTHS
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--obj1')
    parser.add_argument('--obj2')
    parser.add_argument('--base1')
    parser.add_argument('--base2')
    parser.add_argument('--field')

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    obj1 = kwargs['obj1']
    obj2 = kwargs['obj2']

    field1 = kwargs['field']
    field2 = kwargs['field']
    if kwargs['base1'] == 'fico':
        pre1 = 'pkl_efico/' + obj1 + '_' + field1 + '_' + kwargs['base1']
        pre2 = 'pkl_ncorr/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    elif kwargs['base1'] == 'corr':
        pre1 = 'pkl_ecorr/' + obj1 + '_' + field1 + '_' + kwargs['base1']
        pre2 = 'pkl_ncorr/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    elif kwargs['base1'] == 'vary':
        pre1 = 'pkl_evar/' + obj1 + '_' + field1 + '_' + kwargs['base1']
        pre2 = 'pkl_nvary/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    elif kwargs['base1'] == 'fifty':
        pre1 = 'pkl_efifty/' + obj1 + '_' + field1 + '_' + kwargs['base1']
        pre2 = 'pkl_nvary/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    else:
        pre1 = 'pkl_emask/' + obj1 + '_' + field1 + '_' + kwargs['base1']
        pre2 = 'pkl_nmask/' + obj2 + '_' + field2 + '_' + kwargs['base2']

    # pre1 = 'pkls/' + obj1 + '_' + field1 + '_' + kwargs['base1']
    # pre2 = 'nmpkls/' + obj2 + '_' + field2 + '_' + kwargs['base2']

    base1 = '_out.pkl'
    extra1 = pre1 + '_extra' + base1  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res1 = pre1 + '_res' + base1
    sed1 = pre1 + '_sed' + base1
    restwave1 = pre1 + '_restwave' + base1
    spec1 = pre1 + '_spec' + base1
    spswave1 = pre1 + '_spswave' + base1
    chisq1 = pre1 + '_chisq' + base1
    justchi1 = pre1 + '_justchi' + base1

    base2 = '_out.pkl'
    extra2 = pre2 + '_extra' + base2  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res2 = pre2 + '_res' + base2
    sed2 = pre2 + '_sed' + base2
    restwave2 = pre2 + '_restwave' + base2
    spec2 = pre2 + '_spec' + base2
    spswave2 = pre2 + '_spswave' + base2
    chisq2 = pre2 + '_chisq' + base2
    justchi2 = pre2 + '_justchi' + base2

    files1 = [extra1, res1, sed1, restwave1, spec1, spswave1, chisq1, justchi1]
    with open(files1[1], 'rb') as res1:
        results = pickle.load(res1)
    with open(files1[2], 'rb') as sd1:
        sed = pickle.load(sd1)
    with open(files1[3], 'rb') as restwave1:
        wave_rest = pickle.load(restwave1)
    with open(files1[4], 'rb') as spc1:
        spec = pickle.load(spc1)
    with open(files1[5], 'rb') as spswave1:
        sps_wave = pickle.load(spswave1)

    files2 = [extra2, res2, sed2, restwave2, spec2, spswave2, chisq2, justchi2]
    with open(files2[1], 'rb') as res2:
        results2 = pickle.load(res2)
    with open(files2[2], 'rb') as sd2:
        sed2 = pickle.load(sd2)
    with open(files2[3], 'rb') as restwave2:
        wave_rest2 = pickle.load(restwave2)
    with open(files2[4], 'rb') as spc2:
        spec2 = pickle.load(spc2)
    with open(files2[5], 'rb') as spswave2:
        sps_wave2 = pickle.load(spswave2)

    # apply mask
    mask1 = results['obs']['phot_mask']  # mask out -99.0 values
    mask2 = results2['obs']['phot_mask']  # mask out -99.0 values
    fullmask = (mask1 & mask2)
    phot = results['obs']['maggies'][fullmask]
    phot2 = results2['obs']['maggies'][fullmask]
    sed = sed[fullmask]
    sed2 = sed2[fullmask]
    wave_rest = wave_rest[fullmask]
    wave_rest2 = wave_rest2[fullmask]
    print(len(sed), len(sed2))
    spec = spec[110:]
    spec2 = spec2[110:]
    sps_wave = sps_wave[110:]

    # SCALE EELG PHOT SO IT'S MEDIAN = MEDIAN OF SFGs
    med1 = np.median(phot)
    med2 = np.median(phot2)
    f1500_1 = phot[4]
    f1500_2 = phot2[4]
    scale_f1500 = f1500_2 / f1500_1
    scale_phot = med2 / med1
    sedmed = np.median(sed)
    sedmed2 = np.median(sed2)
    scale_sed = sedmed2 / sedmed
    specmed = np.median(spec)
    specmed2 = np.median(spec2)
    scale_spec = specmed2 / specmed
    print(scale_phot, scale_sed, scale_spec)
    phot = scale_phot * phot
    sed *= scale_phot  # scale_sed * sed
    spec *= scale_phot  # 2.3  # scale_spec * spec

    delmod = []
    delspec = []
    delphot = []
    for i in range(len(sed)):
        delmod.append(sed[i] / sed2[i])
    for i in range(len(spec)):
        delspec.append(spec[i] / spec2[i])
    for i in range(len(phot)):
        delphot.append(phot[i] / phot2[i])
    print(len(phot), len(phot2), len(delmod), len(delphot))
    wave_rest = np.asarray(wave_rest)

    # SMOOTH!
    delspec = ndimage.filters.gaussian_filter1d(delspec, sigma=3.)

    # PLOT IT
    # SETUP FIG
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    ax1 = plt.subplot(gs[0])
    # LABEL STUFF
    fs_text = 30
    fs = 20
    fs_ticks = 25
    ax1.set_ylabel(r'Flux Ratio (EELG / SFG)', fontsize=fs_text)
    ax1.set_xlabel(r'Wavelength (Rest) [$\rm \AA$]', fontsize=fs_text)
    # fig.text(0.5, 0.03, r'Wavelength (Rest) [$\rm \AA$]', ha='center', va='bottom', fontsize=fs_text)  # 30

    # GET PLOTTING!
    # ax2 = plt.subplot(gs[1], sharey=ax1, sharex=ax1)
    ymax = 3.25  # 5.75  # 3.9  # 5.75  # 4.1  # 4.5
    ymin = 0.0001  # 0.2
    xmax = 21000
    xmin = 10 ** 3
    ax1.set_ylim(ymin, ymax)  # 3)  # 5)  # (0, 8)
    ax1.set_xlim(xmin, xmax)
    ax1.set_xscale('log')
    # NOT INCLUDING PHOT, MODEL
    # ax1.plot(wave_rest, delmod, marker='D', linestyle='', markerfacecolor='None', markeredgecolor='k',
    #          markersize=10, markeredgewidth=1., label=r'Model')
    # ax1.plot(wave_rest, delphot, marker='o', linestyle='', color='r', label=r'Flux')
    ax1.plot(sps_wave, delspec, color='k', alpha=0.5)  # , label=r'Spectrum')
    ax1.axhline(y=1., linestyle='--', color='k')
    # ax1.legend(numpoints=1, loc='upper left', prop={'size': 20})
    # inset_axes(ax1, width=8 * 0.32, height=8 * 0.28, loc=1)

    # INSET UVJ (EELG)
    ax2 = fig.add_axes([0.741, 0.642, 0.218*0.7, 0.35*0.7])  # left edge, bottom edge, width, height
    objs = [obj1 + '_' + field1, obj2 + '_' + field2]
    uvj.uvj_plot(-1, 'all', objlist=objs, title=False, labels=False, lims=True, size=20, show=False,
                 col=['purple', 'b'], legend=[r'EELG', r'SFG'])
    #     uvj.uvj_plot(-1, 'all', objlist=objs, title=False, labels=False, lims=True, size=25, show=False,

    # EELG, SFG names
    # ax1.text(1.7*10**3, 2.82, str(field1).upper() + '-' + str(obj1) + ', EELG', fontsize=fs_text)
    # ax1.text(1.7*10**3, 2.65, str(field2).upper() + '-' + str(obj2) + ', SFG', fontsize=fs_text)
    # 1.6*10**3 if legend, and legend only if include phot
    ax1.text(1.1*10**3, ymax * 0.94, str(field1).upper() + '-' + str(obj1) + ', EELG', fontsize=fs_text)
    ax1.text(1.1*10**3, ymax * 0.88, str(field2).upper() + '-' + str(obj2) + ', SFG', fontsize=fs_text)

    ax1.axvspan(4800, 5050, color='k', alpha=0.175)  # 0.2

    # TICK PARAMS
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    plt.tick_params(axis='y', which='minor')
    ax1.set_xticks([10 ** 3, 2 * 10 ** 3, 5 * 10 ** 3, 10 ** 4, 2 * 10 ** 4])  # technically works
    ax1.set_xticklabels([r'$10^3$', r'$2\times10^3$', r'$5 \times 10^3$', r'$10^4$', r'$2\times10^4$'], size=fs_ticks)
    ax1.set_yticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])  # technically works
    ax1.set_yticklabels([r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$', r'$2.5$', r'$3.0$'], size=fs_ticks)

    plt.show()
    objs = [obj1, obj2]
    fields = [field1, field2]

'''
Currently running with:
python deltafig.py --obj1=21442 --base1=fico --obj2=7817 --base2=corr --field=cdfs

# python deltafig.py --obj1=21442 --field1=cdfs --base1=fifty --obj2=7817 --field2=cdfs --base2=vary
python deltafig.py --obj1=21442 --base1=fifty --obj2=7817 --base2=vary --field=cdfs
python deltafig.py --obj1=12105 --base1=fifty --obj2=3623 --base2=vary --field=cosmos

# python deltafig.py --obj1=12533 --field1=cdfs --base1=fifty --obj2=7817 --field2=cdfs --base2=vary

python deltafig.py --obj1=15462 --field1=uds --base1=vary --obj2=2920 --field2=cosmos --base2=vary
python deltafig.py --obj1=21076 --field1=cdfs --base1=vary --obj2=7817 --field2=cdfs --base2=vary

'''
