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


def all_plots(fileset, objname, znames, field, font={'fontname': 'Times'}):
    # FIGURES
    fig = plt.figure()

    scale = True
    if scale:
        # FILTERS SCALE FACTOR
        filt_factor1 = 10**8  # 10**7.5  # 7. * 10 ** 7.5
        filt_factor2 = 5 * filt_factor1  # 6. * filt_factor1  # 0.7 * filt_factor1

        # AXIS LIMS
        ymin = 8*10**-2  # 5*10**-2
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
        ymax = 3*10**2  # 2*10**3  # 5*10**3
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
        fs_text = 30
        ylabel = r'Flux [$\mu$Jy]'

    ax1 = plt.subplot(2, 1, 1)  # number cols, number rows, which fig currently on
    # ax1 = plt.subplot(1, 1, 1)  # number cols, number rows, which fig currently on
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlim(xmin, xmax)
    # ax1.set_ylim(10**-0.6, 10**1.3)
    ax1.set_ylim(ymin, ymax)
    # 'axis_name' ('both' is an option), which='both' --> major & minor ticks!
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)

    ax2 = plt.subplot(2, 1, 2, sharex=ax1)
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_xlim(xmin, xmax)
    # ax2.set_ylim(10**-0.6, 10**1.3)
    ax2.set_ylim(ymin2, ymax2)
    ax2.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax2.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)

    for j in range(len(fileset)):
        print(j)
        # PLOT SEPARATE UVJ
        with open(fileset[j][1], 'rb') as res:
            results = pickle.load(res)

        with open(fileset[j][2], 'rb') as sed:
            sed = pickle.load(sed)

        with open(fileset[j][3], 'rb') as restwave:
            wave_rest = pickle.load(restwave)

        with open(fileset[j][4], 'rb') as spec:
            spec = pickle.load(spec)

        with open(fileset[j][5], 'rb') as spswave:
            sps_wave = pickle.load(spswave)

        with open(fileset[j][7], 'rb') as justchi:
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

        zname = znames[j]
        zred = 0
        with open(zname, 'r') as fz:
            hdr_z = fz.readline().split()
        dtype_z = np.dtype([(hdr_z[1], 'S20')] + [(n, np.float) for n in hdr_z[2:]])
        zout = np.loadtxt(zname, comments='#', delimiter=' ', dtype=dtype_z)
        obj_idx = (zout['id'] == objname[j])

        if zout[obj_idx][0][1] >= 0:  # if z_spec exists for this galaxy
            zred = zout[obj_idx][0][1]  # index good for all main cats
        elif objname > 0:
            zred = zout[obj_idx][0][17]  # using zout catalog, z_peak = z_phot; index good for all zout cats

        print(zred)
        if j == 0:
            if scale:
                scale_factor = res_phot[2]  # res_jan[4]  # max(res_jan)
                print(type(res_jan[4]))
                spec_jan = spec_jan / scale_factor
                sed_jan = sed_jan / scale_factor
                res_jan = res_jan / scale_factor
                res_phot = res_phot / scale_factor
                err_jan = err_jan / scale_factor
                err_phot = err_phot / scale_factor
            # ax1.axhline(y=1.)
            # for i in range(len(p_wave)):
            #     if 1425 < p_wave[i] < 1575:
            #         print(i, 'look')
            # ax1.plot(p_wave[2], res_phot[2], '*', markersize=20, color='purple')
            ax1.plot(sps_wave, spec_jan, color='k', alpha=0.5)  # , label=r'Model Spectrum')  # plot spectrum
            ax1.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10, markeredgewidth=1.25,
                     markeredgecolor='k', label=r'Model Photometry')  # plot best fit model
            # ax1.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='purple',
            ax1.errorbar(p_wave, res_phot, yerr=err_phot, marker='o', linestyle='', color='purple',
                         label=r'Observations')  # plot observations
            # ax1.set_ylabel(ylabel, fontsize=fs_text)
            ax1.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            # ax1.axvline(4800, color='k', ls='--')
            # ax1.axvline(5050, color='k', ls='--')
            ax1.axvspan(4800, 5050, color='k', alpha=0.175)  # 0.2
            # ax1.text(700, 3, 'z ~ ' + str(zred) + ', EELG', fontsize=20)
            ax1.text(textx, texty1, str(field[j]).upper() + '-' + str(objname[j]) + ', z = ' + str(zred) + ', EELG',
                     fontsize=fs_text)  # temp setting 100 instead of texty1

            plt.tick_params(axis='y', which='minor')
            # ax1.set_xticks([10**3, 2*10**3, 5*10**3, 10**4, 2*10**4])  # technically works
            # ax1.set_xticklabels([r'$10^3$', r'$2\times10^3$', r'$5 \times 10^3$', r'$10^4$', r'$2\times10^4$'],
            # size=30)
            ax1.set_yticks([10**-2, 10**-1, 10**0, 10**1, 10**2])  # technically works
            ax1.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$'], size=fs_ticks)
            # ax1.set_yticks([10**-0.5, 10**0, 10**0.5, 10**1.])  # technically works
            # ax1.set_yticklabels([r'$10^{-0.5}$', r'$10^0$', r'$10^{0.5}$', r'$10^{1.}$'], size=fs_ticks)

            plt.subplots_adjust(hspace=.0)

            # Redshift for each obj) (note: len(main[0] = 156; elements 153:155 = use, snr, use_nosnr, z_spec)
            widths.fig2(ax1, field[j], zred, scale=(phot.max() * filt_factor1), rest=True)  # WIDTHS
            # widths.plot_filts(ax1, field[j], zred, scale=(phot.max() * filt_factor1 / 2), rest=True)  # WIDTHS

        elif j == 1:
            if scale:
                scale_factor = res_phot[27]  # res_jan[4]  # max(res_jan)
                spec_jan = spec_jan / scale_factor
                sed_jan = sed_jan / scale_factor
                res_jan = res_jan / scale_factor
                res_phot = res_phot / scale_factor
                err_jan = err_jan / scale_factor
                err_phot = err_phot / scale_factor
            # ax2.axhline(y=1.)
            # ax2.plot(p_wave[27], res_phot[27], '*', markersize=20, color='blue')
            ax2.plot(sps_wave, spec_jan, color='k', alpha=0.5)  # , label=r'Model Spectrum')  # plot spectrum
            ax2.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10, markeredgewidth=1.25,
                     markeredgecolor='k', label=r'Model Photometry')  # plot best fit model
            # ax2.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='b',
            ax2.errorbar(p_wave, res_phot, yerr=err_phot, marker='o', linestyle='', color='b',
                         label=r'Observations')  # plot observations
            # ax2.set_ylabel(ylabel, fontsize=fs_text)  # 30
            ax2.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            ax2.axvspan(4800, 5050, color='k', alpha=0.175)  # 0.2
            # ax2.text(700, 3, 'z ~ ' + str(zred) + ', LBG', fontsize=20)
            ax2.text(textx, texty2, str(field[j]).upper() + '-' + str(objname[j]) + ', z = ' + str(zred) + ', SFG',
                     fontsize=fs_text)
            plt.subplots_adjust(hspace=.0)

            # hacking for Vy
            plt.tick_params(axis='y', which='minor')
            ax2.set_xticks([10**3, 2*10**3, 5*10**3, 10**4, 2*10**4])  # technically works
            ax2.set_xticklabels([r'$10^3$', r'$2\times10^3$', r'$5 \times 10^3$', r'$10^4$', r'$2\times10^4$'],
                                size=fs_ticks)
            ax2.set_yticks([10**-2, 10**-1, 10**0, 10**1, 10**2])  # technically works
            ax2.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$'], size=fs_ticks)
            # ax2.set_yticks([10**-0.5, 10**0, 10**0.5, 10**1.])  # technically works
            # ax2.set_yticklabels([r'$10^{-0.5}$', r'$10^0$', r'$10^{0.5}$', r'$10^{1.}$'], size=fs_ticks)
            widths.fig2(ax2, field[j], zred, scale=(phot.max() * filt_factor2), rest=True)  # WIDTHS
            # widths.plot_filts(ax2, field[j], zred, scale=(phot.max() * filt_factor2 / 10**4), rest=True)  # WIDTHS
            plt.subplots_adjust(hspace=.0)
    print('show')

    fig.text(0.07, 0.5, ylabel, fontsize=fs_text, va='center', rotation='vertical', **font)
    plt.xlabel('Wavelength (Rest) [\AA]', fontsize=fs_text, **font)  # 20
    plt.show()

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
    fd = kwargs['field']  # 'cdfs' 'cosmos' 'uds'

    if fd == 'cdfs':
        znames = ['/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout', '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout']
        field1 = 'cdfs'
        field2 = 'cdfs'
    elif fd == 'cosmos':
        znames = ['/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout', '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout']
        field1 = 'cosmos'
        field2 = 'cosmos'
    else:
        znames = ['/home/jonathan/uds/uds.v1.5.8.awk.zout', '/home/jonathan/uds/uds.v1.5.8.awk.zout']
        field1 = 'uds'
        field2 = 'uds'

    if kwargs['base1'] == 'fico':
        pre1 = 'pkl_efico/' + obj1 + '_' + field1 + '_' + kwargs['base1']
        pre2 = 'pkl_ncorr/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    elif kwargs['base1'] == 'corr':
        pre1 = 'pkl_ecorr/' + obj1 + '_' + field1 + '_' + kwargs['base1']
        pre2 = 'pkl_ncorr/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    # pre1 = 'pkl_efico/' + obj1 + '_' + field1 + '_' + kwargs['base1']  # pkl_efico (eelg fifty Myr, corrected)
    # pre2 = 'pkl_ncorr/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    # pre1 = 'pkl_efifty/' + obj1 + '_' + field1 + '_' + kwargs['base1']
    # pre2 = 'pkl_nvary/' + obj2 + '_' + field2 + '_' + kwargs['base2']

    base = '_out.pkl'
    extra1 = pre1 + '_extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res1 = pre1 + '_res' + base
    sed1 = pre1 + '_sed' + base
    restwave1 = pre1 + '_restwave' + base
    spec1 = pre1 + '_spec' + base
    spswave1 = pre1 + '_spswave' + base
    chisq1 = pre1 + '_chisq' + base
    justchi1 = pre1 + '_justchi' + base

    extra2 = pre2 + '_extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res2 = pre2 + '_res' + base
    sed2 = pre2 + '_sed' + base
    restwave2 = pre2 + '_restwave' + base
    spec2 = pre2 + '_spec' + base
    spswave2 = pre2 + '_spswave' + base
    chisq2 = pre2 + '_chisq' + base
    justchi2 = pre2 + '_justchi' + base

    files1 = [extra1, res1, sed1, restwave1, spec1, spswave1, chisq1, justchi1]
    files2 = [extra2, res2, sed2, restwave2, spec2, spswave2, chisq2, justchi2]
    fileset = [files1, files2]
    objs = [obj1, obj2]
    fields = [field1, field2]

    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)

    all_plots(fileset, objs, znames, fields)

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
