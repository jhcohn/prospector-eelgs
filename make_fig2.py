import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import print_sfh
import uvj
import widths  # WIDTHS
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes


def all_plots(fileset, objname, znames, field):
    # FIGURES
    plt.figure()
    # prepares fig to hold two sets side-by-side of: two axes, one atop other, size ratio 3:1

    # FILTERS SCALE FACTOR
    filt_factor = 10 ** 10  # 10 ** 6

    # LEGEND, TEXT LOCATION, FONT SIZE
    loc = 2  # loc=1 (up R), loc=2 (up L), loc=3 (low L), loc=4 (low R); loc=7 (center R)
    textx = 10**4  # 700
    texty = 5*10**2  # 50
    fs = 20
    fs_text = 30

    # AXIS LIMS
    ymin = 2*10**-3  # 10**-7
    ymax = 5*10**3
    xmin = 700  # 600
    xmax = 2.7*10**4  # 3*10**4

    ax1 = plt.subplot(3, 1, 1)  # number cols, number rows, which fig currently on
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    # 'axis_name' ('both' is an option), which='both' --> major & minor ticks!
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)

    ax2 = plt.subplot(3, 1, 2, sharex=ax1)
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax2.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)

    ax3 = plt.subplot(3, 1, 3, sharex=ax1)
    ax3.set_yscale("log")
    ax3.set_xscale("log")
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(ymin, ymax)
    ax3.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax3.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)

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
        phot = results['obs']['maggies'][mask]
        unc = results['obs']['maggies_unc'][mask]
        sed = sed[mask]
        wave_rest = wave_rest[mask]

        # CONVERTING TO microjansky
        factor = 3631 * 10 ** 6
        res_jan = []
        err_jan = []
        sed_jan = []
        for i in range(len(wave_rest)):
            res_jan.append(phot[i] * factor)
            err_jan.append(unc[i] * factor)
            sed_jan.append(sed[i] * factor)
        spec_jan = []
        for i in range(len(spec)):
            spec_jan.append(spec[i] * 3631 * 10 ** 6)

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
        if j == 1:
            ax1.plot(sps_wave, spec_jan, color='k', alpha=0.5)  # , label=r'Model Spectrum')  # plot spectrum
            ax1.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10, markeredgewidth=1.25,
                     markeredgecolor='k', label=r'Model Photometry')  # plot best fit model
            ax1.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='purple',
                         label=r'Observations')  # plot observations
            # ax1.set_ylabel(r'Flux [$\mu$Jy]')
            ax1.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            ax1.axvspan(4800, 5050, color='k', alpha=0.2)
            # ax1.text(700, 3, 'z ~ ' + str(zred) + ', EELG', fontsize=20)
            ax1.text(textx, texty, str(field[j]).upper() + '-' + str(objname[j]) + ', z = ' + str(zred) + ', EELG',
                     fontsize=fs)
            plt.subplots_adjust(hspace=.0)

            # Redshift for each obj) (note: len(main[0] = 156; elements 153:155 = use, snr, use_nosnr, z_spec)
            widths.fig2(ax1, field[j], zred, scale=(phot.max() * filt_factor), rest=True)  # WIDTHS

        elif j == 2:
            ax2.plot(sps_wave, spec_jan, color='k', alpha=0.5)  # , label=r'Model Spectrum')  # plot spectrum
            ax2.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10, markeredgewidth=1.25,
                     markeredgecolor='k', label=r'Model Photometry')  # plot best fit model
            ax2.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='b',
                         label=r'Observations')  # plot observations
            ax2.set_ylabel(r'Flux [$\mu$Jy]', fontsize=fs_text)  # 30
            ax2.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            ax2.axvspan(4800, 5050, color='k', alpha=0.2)
            # ax2.text(700, 3, 'z ~ ' + str(zred) + ', LBG', fontsize=20)
            ax2.text(textx, texty, str(field[j]).upper() + '-' + str(objname[j]) + ', z = ' + str(zred) + ', SFG',
                     fontsize=fs)
            plt.subplots_adjust(hspace=.0)

            widths.fig2(ax2, field[j], zred, scale=(phot.max() * filt_factor), rest=True)  # WIDTHS

        elif j == 0:
            ax3.plot(sps_wave, spec_jan, color='k', alpha=0.5)  # , label=r'Model Spectrum')  # plot spectrum
            ax3.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10, markeredgewidth=1.25,
                     markeredgecolor='k', label=r'Model Photometry')  # plot best fit model
            ax3.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='r',
                         label=r'Observations')  # plot observations
            # ax3.set_ylabel(r'Flux [$\mu$Jy]')
            # ax3.text(700, 1, 'z ~ ' + str(zred) + ', Qui', fontsize=20)
            ax3.text(textx, texty, str(field[j]).upper() + '-' + str(objname[j]) + ', z = ' + str(zred) + ', Qui',
                     fontsize=fs)
            ax3.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            ax3.axvspan(4800, 5050, color='k', alpha=0.2)
            plt.subplots_adjust(hspace=.0)

            widths.fig2(ax3, field[j], zred, scale=(phot.max() * filt_factor), rest=True)  # WIDTHS
    print('show')
    plt.xlabel(r'Wavelength (Rest) [$\AA$]', fontsize=fs_text)  # 20
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--obj1')
    parser.add_argument('--obj2')
    parser.add_argument('--obj3')
    parser.add_argument('--base1')
    parser.add_argument('--base2')
    parser.add_argument('--base3')

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    obj1 = kwargs['obj1']
    obj2 = kwargs['obj2']
    obj3 = kwargs['obj3']

    field1 = 'cdfs'
    field2 = 'cosmos'
    field3 = 'uds'

    pre1 = 'nmpkls/' + obj1 + '_' + field1 + '_' + kwargs['base1']
    # pre2 = 'pkl_emask/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    # pre3 = 'pkl_nmask/' + obj3 + '_' + field3 + '_' + kwargs['base3']
    pre2 = 'pkl_evar/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    pre3 = 'pkl_nvar/' + obj3 + '_' + field3 + '_' + kwargs['base3']
    # pre2 = 'pkls/' + obj2 + '_' + field2 + '_' + kwargs['base2']
    # pre3 = 'nmpkls/' + obj3 + '_' + field3 + '_' + kwargs['base3']

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

    extra3 = pre3 + '_extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res3 = pre3 + '_res' + base
    sed3 = pre3 + '_sed' + base
    restwave3 = pre3 + '_restwave' + base
    spec3 = pre3 + '_spec' + base
    spswave3 = pre3 + '_spswave' + base
    chisq3 = pre3 + '_chisq' + base
    justchi3 = pre3 + '_justchi' + base

    files1 = [extra1, res1, sed1, restwave1, spec1, spswave1, chisq1, justchi1]
    files2 = [extra2, res2, sed2, restwave2, spec2, spswave2, chisq2, justchi2]
    files3 = [extra3, res3, sed3, restwave3, spec3, spswave3, chisq3, justchi3]
    fileset = [files1, files2, files3]
    objs = [obj1, obj2, obj3]
    fields = [field1, field2, field3]

    znames = ['/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout', '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout',
              '/home/jonathan/uds/uds.v1.5.8.awk.zout']

    all_plots(fileset, objs, znames, fields)

'''
Currently running with:
python make_fig2.py --obj1=20752 --obj2=12105 --obj3=5957 --base1=noelg --base2=vary --base3=vary
python make_fig2.py --obj1=20752 --obj2=1824 --obj3=5957 --base1=noelg --base2=fixedmet --base3=noelg
python make_fig2.py --obj1=20752 --obj2=1824 --obj3=5957 --base1=noelg --base2=vary --base3=vary
python make_fig2.py --obj1=20752 --obj2=1824 --obj3=5957 --base1=noelg --base2=newmask --base3=newmask
# try uds 5957 for second obj, since 5206 is an EELG
'''
