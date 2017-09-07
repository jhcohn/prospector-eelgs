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

    ax1 = plt.subplot(3, 1, 1)  # number cols, number rows, which fig currently on
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    # ax1.set_title(field[0] + '-' + objname[0] + ', ' + field[1] + '-' + objname[1] + ', ' + field[2] + '-' +
    # objname[2])

    ax2 = plt.subplot(3, 1, 2, sharex=ax1)
    ax2.set_yscale("log")
    ax2.set_xscale("log")

    ax3 = plt.subplot(3, 1, 3, sharex=ax1)
    ax3.set_yscale("log")
    ax3.set_xscale("log")

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
            ax1.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='r',
                         label=r'Observed Photometry')  # plot observations
            ax1.plot(wave_rest, sed_jan, 'o', color='b', label=r'Model Photometry')  # plot best fit model
            ax1.plot(sps_wave, spec_jan, color='b', alpha=0.5, label=r'Model Spectrum')  # plot spectrum
            ax1.set_ylabel(r'Flux [$\mu$Jy]')
            # ax1.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            ax1.axvspan(4800, 5050, color='k', alpha=0.3)
            ax1.text(700, 3, 'z ~ ' + str(zred) + ', EELG', fontsize=20)
            ax1.text(700, 30, str(field[j]) + '-' + str(objname[j]), fontsize=20)
            plt.subplots_adjust(hspace=.0)

            # Redshift for each obj) (note: len(main[0] = 156; elements 153:155 = use, snr, use_nosnr, z_spec)
            widths.fig2(ax1, field[j], zred, scale=(phot.max() * 10 ** 6), rest=True)  # WIDTHS

        elif j == 2:
            ax2.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='r',
                         label=r'Observed Photometry')  # plot observations
            ax2.plot(wave_rest, sed_jan, 'o', color='b', label=r'Model Photometry')  # plot best fit model
            ax2.plot(sps_wave, spec_jan, color='b', alpha=0.5, label=r'Model Spectrum')  # plot spectrum
            ax2.set_ylabel(r'Flux [$\mu$Jy]')
            # ax2.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            ax2.axvspan(4800, 5050, color='k', alpha=0.3)
            ax2.text(700, 3, 'z ~ ' + str(zred) + ', LBG', fontsize=20)
            ax2.text(700, 30, str(field[j]) + '-' + str(objname[j]), fontsize=20)
            plt.subplots_adjust(hspace=.0)

            widths.fig2(ax2, field[j], zred, scale=(phot.max() * 10 ** 6), rest=True)  # WIDTHS

        elif j == 0:
            ax3.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='r',
                         label=r'Observed Photometry')  # plot observations
            ax3.plot(wave_rest, sed_jan, 'o', color='b', label=r'Model Photometry')  # plot best fit model
            ax3.plot(sps_wave, spec_jan, color='b', alpha=0.5, label=r'Model Spectrum')  # plot spectrum
            ax3.set_ylabel(r'Flux [$\mu$Jy]')
            ax3.text(700, 1, 'z ~ ' + str(zred) + ', Qui', fontsize=20)
            ax3.text(700, 20, str(field[j]) + '-' + str(objname[j]), fontsize=20)
            # ax3.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            ax3.axvspan(4800, 5050, color='k', alpha=0.3)
            plt.subplots_adjust(hspace=.0)

            widths.fig2(ax3, field[j], zred, scale=(phot.max() * 10 ** 6), rest=True)  # WIDTHS
    print('show')
    plt.xlabel(r'Rest frame wavelength [$\AA$]')
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
    pre1 = obj1 + '_' + field1 + '_' + kwargs['base1']
    pre2 = obj2 + '_' + field2 + '_' + kwargs['base2']
    pre3 = obj3 + '_' + field3 + '_' + kwargs['base3']

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
python make_fig2.py --obj1=20752 --obj2=1824 --obj3=5206 --base1=noelg --base2=fixedmet --base3=fixedmet
'''
