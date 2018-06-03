import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import print_sfh
import uvj
import widths  # WIDTHS
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import rc
import stellar_ages as sa
import os


def plot_wrapper(ax, objs, fields, folder, base, color):
    for i in range(len(objs)):
        obj = objs[i]
        fld = fields[i]
        print(obj, i)

        if fld == 'cdfs':
            zname = '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout'  # , '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout']
        elif fld == 'cosmos':
            zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'  # , '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout']
        else:
            zname = '/home/jonathan/uds/uds.v1.5.8.awk.zout'  # , '/home/jonathan/uds/uds.v1.5.8.awk.zout']

        pre = folder + str(obj) + '_' + fld + '_' + base
        fbase = '_out.pkl'
        extra = pre + '_extra' + fbase  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
        res = pre + '_res' + fbase
        sed = pre + '_sed' + fbase
        restwave = pre + '_restwave' + fbase
        spec = pre + '_spec' + fbase
        spswave = pre + '_spswave' + fbase
        chisq = pre + '_chisq' + fbase
        justchi = pre + '_justchi' + fbase
        fileset = [extra, res, sed, restwave, spec, spswave, chisq, justchi]

        all_plots(ax, fileset, obj, zname, color=color)


def all_plots(ax, fileset, objname, zname, color='purple', font={'fontname': 'Times'}):
    print(color)

    home = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
    if os.path.exists(home + fileset[1]):
        with open(fileset[1], 'rb') as res:
            results = pickle.load(res)

        with open(fileset[2], 'rb') as sed:
            sed = pickle.load(sed)

        with open(fileset[3], 'rb') as restwave:
            wave_rest = pickle.load(restwave)

        with open(fileset[4], 'rb') as spec:
            spec = pickle.load(spec)

        with open(fileset[5], 'rb') as spswave:
            sps_wave = pickle.load(spswave)

        with open(fileset[7], 'rb') as justchi:
            chi = pickle.load(justchi)

        # MASK
        # create mask
        mask = results['obs']['phot_mask']  # mask out -99.0 values
        wave_rest = np.asarray(wave_rest)
        # no longer need this once I use new output that creates wave_rest as an array
        # ly_mask = (1180 < wave_rest) & (wave_rest < 1260)  # mask out ly-alpha values
        # mask[ly_mask] = False  # combine the masks!

        # apply mask
        phot = results['obs']['maggies']# [mask]
        unc = results['obs']['maggies_unc']# [mask]
        p_wave = wave_rest
        sed = sed[mask]
        wave_rest = wave_rest[mask]
        print(len(phot), len(p_wave), 'meeee')

        # CONVERTING TO microjansky
        factor = 3631 * 10 ** 6
        res_jan = []
        err_jan = []
        sed_jan = []
        res_phot = []
        err_phot = []
        for i in range(len(wave_rest)):
            res_jan.append(phot[i] * factor)
            err_jan.append(unc[i] * factor)
            sed_jan.append(sed[i] * factor)
        spec_jan = []
        for i in range(len(spec)):
            spec_jan.append(spec[i] * 3631 * 10 ** 6)
        for i in range(len(p_wave)):
            res_phot.append(phot[i] * factor)
            err_phot.append(unc[i] * factor)

        zred = 0
        with open(zname, 'r') as fz:
            hdr_z = fz.readline().split()
        dtype_z = np.dtype([(hdr_z[1], 'S20')] + [(n, np.float) for n in hdr_z[2:]])
        zout = np.loadtxt(zname, comments='#', delimiter=' ', dtype=dtype_z)
        obj_idx = (zout['id'] == str(objname))

        spec_z = False
        if zout[obj_idx][0][1] >= 0:  # if z_spec exists for this galaxy
            zred = zout[obj_idx][0][1]  # index good for all main cats
            spec_z = True
        elif objname > 0:
            zred = zout[obj_idx][0][17]  # using zout catalog, z_peak = z_phot; index good for all zout cats

        print(zred)
        '''
        if scale:
            scale_factor = res_phot[2]  # res_jan[4]  # max(res_jan)
            print(type(res_jan[4]))
            spec_jan = spec_jan / scale_factor
            sed_jan = sed_jan / scale_factor
            res_jan = res_jan / scale_factor
            res_phot = res_phot / scale_factor
            err_jan = err_jan / scale_factor
            err_phot = err_phot / scale_factor
        '''
        # ax.plot(sps_wave, spec_jan, color='k', alpha=0.5)  # , label=r'Model Spectrum')  # plot spectrum
        spec_jan /= np.median(spec_jan)
        ax.plot(sps_wave, spec_jan, color=color, alpha=0.5)  # , label=r'Model Spectrum')  # plot spectrum
        # ax.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10, markeredgewidth=1.25,
        #         markeredgecolor='k')  # , label=r'Model Photometry')  # plot best fit model
        # ax.errorbar(p_wave, res_phot, yerr=err_phot, marker='o', linestyle='', color=color)  # ,
                    # label=r'Observations')  # plot observations
        # ax.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
        '''
        if spec_z:
            ax1.text(textx, texty1, str(field[j]).upper() + '-' + str(objname[j]) + r', z$_{\rm spec}$ = '
                     + str(zred) + ', EELG', fontsize=fs_text)
        else:
            ax1.text(textx, texty1, str(field[j]).upper() + '-' + str(objname[j]) + r', z$_{\rm phot}$ = ' +
                     str(round(zred, 2)) + ', EELG', fontsize=fs_text)
        '''
    else:
        print('NOOOO', fileset[1])


if __name__ == "__main__":
    fico = True
    if fico:
        folders = ['pkl_efico/', 'pkl_nfico/']
        base = 'fico'

    # FIGURES
    fig = plt.figure()

    scale = False  # True
    if scale:
        # FILTERS SCALE FACTOR
        filt_factor1 = 10**8  # 10**7.5  # 7. * 10 ** 7.5
        filt_factor2 = 5 * filt_factor1  # 6. * filt_factor1  # 0.7 * filt_factor1

        # AXIS LIMS
        ymin = 2*10**-2  # 5*10**-2
        ymax = 400  # 300
        xmin = 700  # 600
        xmax = 2.7 * 10 ** 4  # 3*10**4
        ymin2 = ymin
        ymax2 = ymax

        # LEGEND, TEXT LOCATION, FONT SIZE
        loc = 2  # loc=1 (up R), loc=2 (up L), loc=3 (low L), loc=4 (low R); loc=7 (center R)
        textx = 720  # 1.5 * 10 ** 3  # 10**4  # 700
        texty1 = 10  # 30  # 20
        texty2 = texty1
        fs = 20
        fs_ticks = 25
        fs_text = 30
        ylabel = r'Scaled F$_\nu$ [$1500$ \AA]'

    else:
        # FILTERS SCALE FACTOR
        filt_factor1 = 10 ** 8.5  # 10 ** 10  # 10 ** 6
        filt_factor2 = 1.75 * 10 ** 8.5  # 10 ** 10  # 10 ** 6

        # AXIS LIMS
        ymin = 3*10**-2  # 2*10**-3  # 10**-7
        ymax = 7*10**2  # 3*10**2  # 2*10**3  # 5*10**3
        xmin = 700  # 600
        xmax = 2.7*10**4  # 3*10**4
        ymin2 = 9*10**-2  # ymin
        ymax2 = 9*10**2  # ymax

        # LEGEND, TEXT LOCATION, FONT SIZE
        loc = 2  # loc=1 (up R), loc=2 (up L), loc=3 (low L), loc=4 (low R); loc=7 (center R)
        textx = 1.7*10**3  # 10**4  # 700
        texty1 = 10**2  # 50
        texty2 = 3*10**2
        fs = 20
        fs_ticks = 25
        fs_text = 30
        ylabel = r'Flux [$\mu$Jy]'

    e_objs, e_fields, l_objs, l_fields = sa.get_gal_lists(base='fico', objlists=True)

    ax1 = plt.subplot(1, 2, 1)  # number cols, number rows, which fig currently on
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs) # 'axis_name', which='both' --> major & minor!
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    ax1.tick_params(axis='x', which='major', pad=10)
    ax1.tick_params(axis='y', which='minor')
    ax1.set_xticks([10 ** 3, 2 * 10 ** 3, 5 * 10 ** 3, 10 ** 4, 2 * 10 ** 4])  # technically works
    ax1.axvspan(4800, 5050, color='k', alpha=0.175)  # 0.2
    ax1.set_xticklabels([r'$10^3$', r'$2\times10^3$', r'$5 \times 10^3$', r'$10^4$', r'$2\times10^4$'],
                        size=fs_ticks)
    ax1.set_yticks([10 ** -1, 10 ** 0, 10 ** 1, 10 ** 2])  # technically works  # 10**-2,
    ax1.set_yticklabels([r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$'], size=fs_ticks)  # r'$10^{-2}$',
    plot_wrapper(ax1, e_objs, e_fields, folders[0], base, color='purple')

    ax2 = plt.subplot(1, 2, 2, sharey=ax1)
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax2.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    plt.subplots_adjust(hspace=.0)
    ax2.tick_params(axis='y', which='minor')
    ax2.tick_params(axis='x', which='major', pad=10)
    ax2.set_xticks([10 ** 3, 2 * 10 ** 3, 5 * 10 ** 3, 10 ** 4, 2 * 10 ** 4])  # technically works
    ax2.set_xticklabels([r'$10^3$', r'$2\times10^3$', r'$5 \times 10^3$', r'$10^4$', r'$2\times10^4$'],
                        size=fs_ticks)
    # ax2.set_yticks([10 ** -1, 10 ** 0, 10 ** 1, 10 ** 2])  # technically works  # 10**-2,
    # ax2.set_yticklabels([])  # [r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$'], size=fs_ticks)  # r'$10^{-2}$',
    ax2.axvspan(4800, 5050, color='k', alpha=0.175)  # 0.2
    plot_wrapper(ax2, l_objs, l_fields, folders[1], base, color='b')

    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    font = {'fontname': 'Times'}

    plt.subplots_adjust(wspace=.0)
    plt.subplots_adjust(hspace=.0)

    fig.text(0.07, 0.5, ylabel, fontsize=fs_text, va='center', rotation='vertical', **font)
    fig.text(0.5, 0.02, 'Wavelength (Rest) [\AA]', ha='center', fontsize=fs_text, **font)  # 20
    # plt.xlabel('Wavelength (Rest) [\AA]', fontsize=fs_text, **font)  # 20
    plt.show()

'''
Currently running with:
python new_fig1.py --obj1=21442 --base1=fico --obj2=7817 --base2=corr --field=cdfs

python new_fig1.py --obj1=21442 --base1=fifty --obj2=7817 --base2=vary --field=cdfs
python new_fig1.py --obj1=12105 --base1=fifty --obj2=2729 --base2=vary --field=cosmos  # pretty good, deltafig 2xtreme
python new_fig1.py --obj1=12105 --base1=fifty --obj2=3623 --base2=vary --field=cosmos  # pretty good, deltafig 2xtreme!!
python new_fig1.py --obj1=12105 --base1=fifty --obj2=8233 --base2=vary --field=cosmos  # not bad...
python new_fig1.py --obj1=12105 --base1=fifty --obj2=12535 --base2=vary --field=cosmos  # pretty good, big-ish z-diff
python new_fig1.py --obj1=12105 --base1=fifty --obj2=12676 --base2=vary --field=cosmos  # not bad, deltafig good emlines
# but not great slope

# WAS DOING:
python new_fig1.py --obj1=12533 --base1=fifty --obj2=7817 --base2=vary
# look for smaller photometric errors than 12533
# Shape-wise, I like 11462. 12903 not as bad error-wise. 15124 best errors, but lack of impressive emission points.
# 18561 not super terrible. 21442 quite good.

Prior to taking this ratio, I scaled the EELG such that its median flux is the same as the median flux of the SFG (the
SFG inherently has a larger continuum, the only purpose of this scaling is so that a ratio of "1" corresponds to
effectively equivalent fluxes on the EELG and SFG compared to their respective continua; the scaling doesn't change the
actual shape of the plot at all).
'''
