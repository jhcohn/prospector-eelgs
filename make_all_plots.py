import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import print_sfh
import uvj
import widths  # WIDTHS
from matplotlib import gridspec


def all_plots(files, objname, field, kwbase, loc='upper left', sep_uvj=False, curves=False, ujy=True):
    # PLOT SEPARATE UVJ
    if sep_uvj:
        uvj.uvj_plot(objname, field)

    # PLOT SFH, SED+
    # files should be a list of output files in order: sfh, res, sed, restwave, spec, spswave, chisq, justchi
    print_sfh.plotter(files[0])
    print_sfh.plotter(files[0], specific=True)

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

    # FIGURES
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])  # prepares fig to hold two axes, one atop other, size ratio 3:1

    ax1 = plt.subplot(gs[0])  # ax1 = bigger upper axis
    # plt.axvline(x=5270, color='k')  # proving no offset in two axes
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_title(field + '-' + obj)

    # MASK
    # create mask
    mask = results['obs']['phot_mask']  # mask out -99.0 values
    wave_rest = np.asarray(wave_rest)  # no longer need this once I use new output that creates wave_rest as an array
    ly_mask = (1180 < wave_rest) & (wave_rest < 1260)  # mask out ly-alpha values
    mask[ly_mask] = False  # combine the masks!

    # apply mask
    phot = results['obs']['maggies'][mask]
    unc = results['obs']['maggies_unc'][mask]
    sed = sed[mask]
    wave_rest = wave_rest[mask]
    chi = chi[mask]

    if ujy:
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

        ax1.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='r',
                     label=r'Observed Photometry')  # plot observations
        ax1.plot(wave_rest, sed_jan, 'o', color='b', label=r'Model Photometry')  # plot best fit model
        ax1.plot(sps_wave, spec_jan, color='b', alpha=0.5, label=r'Model Spectrum')  # plot spectrum
        ax1.set_ylabel(r'$\mu$Jy')
    else:
        # MAGGIES
        ax1.errorbar(wave_rest, phot, yerr=unc, marker='o', linestyle='', color='r', label=r'Observed Photometry')
        # plot observations^
        ax1.plot(wave_rest, sed, 'o', color='b', label=r'Model Photometry')  # plot best fit model
        ax1.plot(sps_wave, spec, color='b', alpha=0.5, label=r'Model Spectrum')  # plot spectrum
        ax1.set_ylabel(r'Maggies')

    # OVERLAY SED WITH FILTER WIDTHS
    if curves:
        # Redshift for each obj) (note: len(main[0] = 156; elements 153:155 = use, snr, use_nosnr, z_spec)
        zred = 0
        if field == 'cdfs':
            zname = '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout'
        elif field == 'cosmos':
            zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'  # redshift catalog
        elif field == 'uds':
            zname = '/home/jonathan/uds/uds.v1.5.8.awk.zout'

        with open(zname, 'r') as fz:
            hdr_z = fz.readline().split()
        dtype_z = np.dtype([(hdr_z[1], 'S20')] + [(n, np.float) for n in hdr_z[2:]])
        zout = np.loadtxt(zname, comments='#', delimiter=' ', dtype=dtype_z)
        obj_idx = (zout['id'] == objname)

        if zout[obj_idx][0][1] >= 0:  # if z_spec exists for this galaxy
            zred = zout[obj_idx][0][1]  # index good for all main cats
        elif objname > 0:
            zred = zout[obj_idx][0][17]  # using zout catalog, z_peak = z_phot; index good for all zout cats

        widths.some_filts(field, zred, scale=(phot.max() / 10 ** 2), rest=True)  # WIDTHS

    ax2 = plt.subplot(gs[1], sharex=ax1)  # ax2 = smaller lower axis
    line2 = ax2.plot(wave_rest, chi, 'o', color='k')  # plot chi
    plt.axhline(y=0, color='k')  # plot horizontal line at y=0 on chi plot
    plt.setp(ax1.get_xticklabels(), visible=False)  # hide xtick labels on upper axis
    yticks = ax2.yaxis.get_major_ticks()  # show ytick labels on lower axis
    yticks[-1].label1.set_visible(False)  # hide uppermost ytick label on lower axis to prevent overlap
    ax2.set_ylabel(r'$\chi$')
    ax2.set_xlabel(r'Rest frame wavelength')

    ax1.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
    # (line11, line12, line13), labels=[r'Model', r'Spectrum', r'Observed Photometry'],

    # plt.axvline(x=5270, color='k')  # proving no offset in two axes
    plt.subplots_adjust(hspace=.0)

    if not sep_uvj and not curves:  # as long as not plotting a separate uvj plot
        # UVJ inset on SED+ plot
        from mpl_toolkits.axes_grid.inset_locator import inset_axes
        loc_uvj = 1
        if kwbase == 'noelg':
            loc_uvj = 7
        inset_axes(ax1, width="20%", height=2., loc=loc_uvj)  # create inset axis: width (%), height (inches), location
        # loc=1 (upper right), loc=2 (upper left) --> loc=3 (lower left), loc=4 (lower right); loc=7 (center right)
        # https://stackoverflow.com/questions/10824156/matplotlib-legend-location-numbers
        uvj.uvj_plot(objname, field, title=False, labels=False, lims=True, size=20)  # add uvj plot to inset axis

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--obj')
    parser.add_argument('--field')
    parser.add_argument('--base')

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    obj = kwargs['obj']

    field = kwargs['field']
    pre = obj + '_' + field + '_' + kwargs['base']

    base = '_out.pkl'
    extra = pre + '_extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res = pre + '_res' + base
    sed = pre + '_sed' + base
    restwave = pre + '_restwave' + base
    spec = pre + '_spec' + base
    spswave = pre + '_spswave' + base
    chisq = pre + '_chisq' + base
    justchi = pre + '_justchi' + base

    files = [extra, res, sed, restwave, spec, spswave, chisq, justchi]

    all_plots(files, obj, field, kwargs['base'], loc='upper left', curves=False)

'''
Currently running with:
python make_all_plots.py --obj=7730 --field=cosmos --base=fixedmet
'''
