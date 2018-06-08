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
import pandas
from scipy import ndimage


def distribution_wrapper(ax, objs, fields, folder, base, color, label):
    bounds = [750, 800, 850, 900, 950, 10**3, 1.1*10**3, 1.2*10**3, 1.3*10**3, 1.4*10**3, 1.5*10**3, 1.6*10**3,
              1.7*10**3, 1.8*10**3, 1.9*10**3, 2*10**3, 2.2*10**3, 2.4*10**3, 2.6*10**3, 2.8*10**3, 3.*10**3,
              3.3*10**3, 3.6*10**3, 3.9*10**3, 4.2*10**3, 4.5*10**3, 4.8*10**3, 5.1*10**3, 5.5*10**3,
              7.5*10**3, 8.*10**3, 9.*10**3, 10**4, 1.1*10**4]
              # 6.*10**3, 6.5*10**3, 7.*10**3,
              # , 1.2*10**4, 1.4*10**4,
              #1.6*10**4, 1.8*10**4, 2.*10**4, 2.25*10**4, 2.5*10**4]
    waves = []
    for x in range(len(bounds) - 1):
        waves.append((bounds[x] + bounds[x+1]) / 2)
    totes = np.zeros(shape=(len(bounds) - 1, len(objs)*40))

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

        wave1, phot1 = all_plots(ax, fileset, obj, zname, color=color, do_plot=False)

        # NEW: plot average values across wavelengths
        for wv in range(len(wave1)):
            for bn in range(len(bounds) - 1):
                if bounds[bn] <= wave1[wv] <= bounds[bn + 1]:
                    totes[bn, i] = phot1[wv]

    plus = []
    avgs = []
    minus = []
    for to in range(len(totes)):  # at each wave point (within each boundary)
        summer = 0.
        typicals = []
        for tes in range(len(totes[to])):  # for each point that fell within that boundary
            if totes[to, tes] != 0.:
                summer += 1
                typicals.append(totes[to, tes])
        if summer >= 3:
            avgs.append(np.percentile(typicals, [16, 50., 84.])[1])
            plus.append(np.percentile(typicals, [16, 50., 84.])[2])
            minus.append(np.percentile(typicals, [16, 50., 84.])[0])
        else:
            avgs.append(0.)
            plus.append(0.)
            minus.append(0.)

    #print(len(avgs), len(plus), len(minus))
    print(avgs)
    print(len(waves), waves)
    ind_to_pop = []
    spec_pop = []
    for pl in range(len(avgs)):
        if minus[pl] <= 0.01:
            ind_to_pop.append(pl)
            # if avgs[pl] <= 0.05 and avgs[pl] != 0.:
            #     spec_pop.append(pl)
        # elif minus[pl] <= 3*10**-2:
        #     ind_to_pop.append(pl)

    # print(ind_to_pop)
    # print([waves[x] for x in spec_pop])
    # waves_sp = [waves[x] for x in spec_pop]
    # up_lims = [plus[x] for x in spec_pop]

    for ind in sorted(ind_to_pop, reverse=True):
        avgs.pop(ind)
        plus.pop(ind)
        minus.pop(ind)
        waves.pop(ind)
    # print(waves, avgs)

    split_waves1 = []
    avgs1, avgs2, avgsm, minus1, minus2, minusm, plus1, plus2, plusm = [], [], [], [], [], [], [], [], []
    split_waves2 = []
    mid_waves = []
    for wav in range(len(waves)):
        if waves[wav] <= 4.5*10**3:
            split_waves1.append(waves[wav])
            avgs1.append(avgs[wav])
            minus1.append(minus[wav])
            plus1.append(plus[wav])
        elif waves[wav] >= 5.5*10**3:
            split_waves2.append(wav)
            avgs2.append(avgs[wav])
            minus2.append(minus[wav])
            plus2.append(plus[wav])
        else:
            mid_waves.append(wav)
            avgsm.append(avgs[wav])
            minusm.append(minus[wav])
            plusm.append(plus[wav])

    # '''
    avgs = pandas.Series((x for x in avgs))
    minus = pandas.Series((x for x in minus))
    plus = pandas.Series((x for x in plus))
    avgs = pandas.rolling_median(avgs, 2)
    minus = pandas.rolling_median(minus, 2)
    plus = pandas.rolling_median(plus, 2)
    '''
    avgs1 = pandas.Series((x for x in avgs1))
    minus1 = pandas.Series((x for x in minus1))
    plus1 = pandas.Series((x for x in plus1))
    avgs1 = pandas.rolling_median(avgs1, 2)
    minus1 = pandas.rolling_median(minus1, 2)
    plus1 = pandas.rolling_median(plus1, 2)
    avgs2 = pandas.Series((x for x in avgs2))
    minus2 = pandas.Series((x for x in minus2))
    plus2 = pandas.Series((x for x in plus2))
    avgs2 = pandas.rolling_median(avgs2, 2)
    minus2 = pandas.rolling_median(minus2, 2)
    plus2 = pandas.rolling_median(plus2, 2)
    avgs = avgs1.tolist() + avgsm + avgs2.tolist()
    plus = plus1.tolist() + plusm + plus2.tolist()
    minus = minus1.tolist() + minusm + minus2.tolist()
    '''

    print(len(waves), waves)
    ax.fill_between(waves, minus, plus, color=color, alpha=0.3)  # fill region between +/- 1sigma
    ax.plot(waves, avgs, '-', lw=2, color=color, label=label)  # , markersize=20, label=label)
    ax.plot(waves, plus, '-', lw=2, color=color)  # , markersize=20, label=label)
    ax.plot(waves, minus, '-', lw=2, color=color)  # , markersize=20, label=label)
    # ax.plot(waves_sp, up_lims, 'v', color=color, markersize=20)


def plot_wrapper(ax, objs, fields, folder, base, color, label):
    bounds = [700, 750, 800, 850, 900, 950, 10**3, 1.1*10**3, 1.2*10**3, 1.3*10**3, 1.4*10**3, 1.5*10**3, 1.6*10**3,
              1.7*10**3, 1.8*10**3, 1.9*10**3, 2*10**3, 2.2*10**3, 2.4*10**3, 2.6*10**3, 2.8*10**3, 3.*10**3,
              3.3*10**3, 3.6*10**3, 3.9*10**3, 4.2*10**3, 4.5*10**3, 4.8*10**3, 5.1*10**3, 5.5*10**3,
              6.*10**3, 6.5*10**3, 7.*10**3, 7.5*10**3, 8.*10**3, 9.*10**3, 10**4, 1.1*10**4, 1.2*10**4, 1.4*10**4,
              1.6*10**4, 1.8*10**4, 2.*10**4, 2.2*10**4, 2.4*10**4, 2.7*10**4]
    waves = []
    for x in range(len(bounds) - 1):
        waves.append((bounds[x] + bounds[x+1]) / 2)
    totes = np.zeros(shape=(len(bounds) - 1, len(objs)*40))

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

        wave1, phot1 = all_plots(ax, fileset, obj, zname, color=color)

        # NEW: plot average values across wavelengths
        for wv in range(len(wave1)):
            for bn in range(len(bounds) - 1):
                if bounds[bn] <= wave1[wv] <= bounds[bn + 1]:
                    totes[bn, i] = phot1[wv]

    avgs = []
    for to in range(len(totes)):  # at each wave point (within each boundary)
        summer = 0.
        typicals = []
        for tes in range(len(totes[to])):  # for each point that fell within that boundary
            if totes[to, tes] != 0.:
                summer += 1
                typicals.append(totes[to, tes])
        if summer >= 4:
            avgs.append(np.percentile(typicals, [16, 50., 84.])[1])
        else:
            avgs.append(0.)

    print(waves, avgs)
    ax.plot(waves, avgs, '*', color=color, markersize=20, label=label)


def all_plots(ax, fileset, objname, zname, color='purple', do_plot=True, font={'fontname': 'Times'}):
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
        #spec_jan /= np.median(spec_jan)
        #ax.plot(sps_wave, spec_jan, color=color, alpha=0.5)  # , label=r'Model Spectrum')  # plot spectrum
        res_phot /= np.median(res_phot)

        res_phot = [x for _,x in sorted(zip(p_wave, res_phot))]
        p_wave = sorted(p_wave)
        if color == 'purple':
            alph = 0.2
        else:
            alph = 0.07
        if do_plot:
            ax.plot(p_wave, res_phot, 'o', color=color, alpha=alph)
        '''
        # ax.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10, markeredgewidth=1.25,
        #         markeredgecolor='k')  # , label=r'Model Photometry')  # plot best fit model
        # ax.errorbar(p_wave, res_phot, yerr=err_phot, marker='o', linestyle='', color=color)  # ,
                    # label=r'Observations')  # plot observations
        # ax.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
        if spec_z:
            ax1.text(textx, texty1, str(field[j]).upper() + '-' + str(objname[j]) + r', z$_{\rm spec}$ = '
                     + str(zred) + ', EELG', fontsize=fs_text)
        else:
            ax1.text(textx, texty1, str(field[j]).upper() + '-' + str(objname[j]) + r', z$_{\rm phot}$ = ' +
                     str(round(zred, 2)) + ', EELG', fontsize=fs_text)
        '''
        return p_wave, res_phot
    else:
        print('NOOOO', fileset[1])
        return [0.], [0.]


if __name__ == "__main__":
    fico = 1
    comp_fast = 0

    distrib = 1
    log = 0

    if fico:
        folders = ['pkl_efico/', 'pkl_nfico/']
        base1 = 'fico'
        base2 = 'fico'
    elif comp_fast:
        folders = ['pkl_efico/', 'pkl_efast/']
        base1 = 'fico'
        base2 = 'fast'

    # FIGURES
    fig = plt.figure()
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    font = {'fontname': 'Times'}

    # LEGEND, TEXT LOCATION, FONT SIZE
    loc = 2  # loc=1 (up R), loc=2 (up L), loc=3 (low L), loc=4 (low R); loc=7 (center R)
    textx = 1.7*10**3  # 10**4  # 700
    texty1 = 10**2  # 50
    texty2 = 3*10**2
    fs = 25  # 20
    fs_ticks = 30  # 25
    fs_text = 35  # 30
    ylabel = r'F$_\nu$ [scaled]'  # [$\mu$Jy]'

    e_objs, e_fields, l_objs, l_fields = sa.get_gal_lists(base=base1, objlists=True)

    ax1 = plt.subplot(1, 1, 1)  # number cols, number rows, which fig currently on
    if log:
        ax1.set_yscale("log")
        ax1.set_xscale("log")
        # AXIS LIMS
        ymin = 9 * 10 ** -2  # 3*10**-2  # 2*10**-3  # 10**-7
        ymax = 1.5 * 10 ** 1  # 7*10**2  # 3*10**2  # 2*10**3  # 5*10**3
        xmin = 800  # 700  # 600
        xmax = 2.5 * 10 ** 4  # 3*10**4
    else:
        # AXIS LIMS
        ax1.set_xscale("log")
        ymin = 0.  # 2*10**-3  # 10**-7
        ymax = 7.  # 15  # 7*10**2  # 3*10**2  # 2*10**3  # 5*10**3  # NOTE: use 15 if not doing percs
        xmin = 10**3  # 800  # 600
        xmax = 10**4  # 2.5 * 10 ** 4  # 3*10**4
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs) # 'axis_name', which='both' --> major & minor!
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    ax1.tick_params(axis='x', which='major', pad=10)
    ax1.tick_params(axis='y', which='minor')
    if log:
        ax1.set_yticks([10 ** -1, 10 ** 0, 10 ** 1])  # , 10 ** 2])  # technically works  # 10**-2,
        ax1.set_yticklabels([r'$10^{-1}$', r'$10^0$', r'$10^1$'], size=fs_ticks)  # r'$10^{-2}$',  , r'$10^2$'
    else:
        # ax1.set_xticks([2000, 4000, 6000, 8000, 10**4, 1.2*10**4, 1.4*10**4, 1.6*10**4, 1.8*10**4, 2*10**4, 2.2*10**4,
        #                 2.4*10**4])  # technically works
        # ax1.set_xticklabels([r'$2\times10^3$', r'$4\times10^3$', r'$6\times10^3$', r'$8\times10^3$', r'$10^4$',
        #                      r'$1.2\times10^4$', r'$1.4\times10^4$', r'$1.6\times10^4$', r'$1.8\times10^4$',
        #                      r'$2\times10^4$', r'$2.2\times10^4$', r'$2.4\times10^4$'], size=fs_ticks)  # , r'$2.6\times10^4$'
        ax1.set_yticks([0., 2., 4., 6.])
        ax1.set_yticklabels([r'$0$', r'$2$', r'$4$', r'$6$'], size=fs_ticks)  # , r'$8$', r'$10$', r'$12$', r'$14$'

    ax1.set_xticks([10 ** 3, 2 * 10 ** 3, 5 * 10 ** 3, 10 ** 4])  # technically works  #  2 * 10 ** 4
    ax1.axvspan(4800, 5050, color='k', alpha=0.25)  # 0.2
    ax1.set_xticklabels([r'$10^3$', r'$2\times10^3$', r'$5 \times 10^3$', r'$10^4$'], size=fs_ticks)  # , r'$2\times10^4$'

    if comp_fast:
        plot_wrapper(ax1, e_objs, e_fields, folders[1], base2, color='b', label='FAST-like EELGs')
        plot_wrapper(ax1, e_objs, e_fields, folders[0], base1, color='purple', label='Regular EELGs')
        fig.text(0.2, 0.8, 'Regular EELGs', fontsize=fs_text, **font)
        fig.text(0.6, 0.8, 'FAST-like runs', fontsize=fs_text, **font)
    elif distrib:
        distribution_wrapper(ax1, l_objs, l_fields, folders[1], base2, color='b', label='SFGs')
        distribution_wrapper(ax1, e_objs, e_fields, folders[0], base1, color='purple', label='EELGs')
        '''
        andromeda = [7.65, 7.95, 8.75, 9.3, 9.55, 9.65, 9.7,
                     9., 8.65, 8.5, 8.2]
        an_waves = [1800, 2200, 3400, 4800, 6000, 7500, 9000,
                    35000, 45000, 59000, 80000]
        andromeda /= np.median(andromeda)
        ax1.plot(an_waves, andromeda, '--', lw=2, color='k', label=r'Andromeda (Groves+2012)')
        '''
    else:
        plot_wrapper(ax1, l_objs, l_fields, folders[1], base2, color='b', label='SFGs')
        plot_wrapper(ax1, e_objs, e_fields, folders[0], base1, color='purple', label='EELGs')

    ax1.legend(numpoints=1, loc=loc, prop={'size': fs_text})  # , line2) ... , r'$\chi$']

    # plt.subplots_adjust(wspace=.0)
    # plt.subplots_adjust(hspace=.0)

    ax1.tick_params(axis='x', which='major', pad=15)
    fig.text(0.575, 0.82, r'H$\beta$+[OIII]', fontsize=fs_text, **font)  # 0.575, 0.85
    fig.text(0.04, 0.5, ylabel, fontsize=fs_text, va='center', rotation='vertical', **font)  # 0.08
    # fig.text(0.5, 0.92, 'Wavelength (Rest) [\AA]', ha='center', fontsize=fs_text, **font)  # 0.02
    plt.xlabel('Wavelength (Rest) [\AA]', fontsize=fs_text, **font)  # 20
    plt.show()
