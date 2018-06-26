import os
import make_all_plots as map
import pickle
import numpy as np
import print_sfh
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import uvj

home = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
eelgs = 0  # True=1
tt = 0
tfast = 1
sfh = 0  # if not sfh, then seds
chi_stuff = 0
folders = ['pkl_tt', 'pkl_tfn']  # ['pkl_masstest', 'pkl_ncorr']
# ['pkl_efix', 'pkl_nvary']  # ['pkl_efifty', 'pkl_nvary']  # ['pkl_evar', 'pkl_nvary']  # ['pkl_emask', 'pkl_nmask']
file = 'tfn'  # 'tt'  # 'masstest'  # 'fifty'  # 'vary'  # 'newmask'
fs = 20  # 30
fs_text = 30
textx = 1100
texty = 30  # 300  # 500

objs = []
fields = []
count = 0
if eelgs:
    folder = folders[0]
    galxs = open(home + 'eelg_specz_ids1', 'r')  # open(base + 'specz_ids', 'r')
    for line in galxs:
        skip = False
        if line[0] == '#':  # ignore commented lines
            skip = True
        space = 0  # counts which column we're in, using python indexing (space = 0 --> 1st column)
        cols = line.split()
        # field = cols[0]
        if not skip:
            if int(cols[0]) - 200000 > 0:
                field = 'uds'
                obj = str(int(cols[0]) - 200000)
            elif int(cols[0]) - 100000 > 0:
                field = 'cosmos'
                obj = str(int(cols[0]) - 100000)
            else:
                field = 'cdfs'
                obj = cols[0]
            count += 1
            objs.append(obj)
            fields.append(field)
    galxs.close()
    # print(objs)
elif tt:
    # galxs = ['5519_cdfs_tt', '5593_cdfs_tt', '547_cdfs_tt']
    folder = 'pkl_tt'
    objs = ['5519', '5593', '5475']
    fields = ['cdfs', 'cdfs', 'cdfs']
elif tfast:
    # galxs = ['5519_cdfs_tt', '5593_cdfs_tt', '547_cdfs_tt']
    folder = 'pkl_tfn'
    objs = ['5519', '5593', '5475']
    fields = ['cdfs', 'cdfs', 'cdfs']
else:
    folder = folders[1]
    galxs = open(home + 'lbg_ids1', 'r')
    for line in galxs:
        count += 1
        if int(line) - 200000 > 0:
            objs.append(str(int(line) - 200000))
            fields.append('uds')
        elif int(line) - 100000 > 0:
            objs.append(str(int(line) - 100000))
            fields.append('cosmos')
        else:
            objs.append(str(int(line)))
            fields.append('cdfs')
    galxs.close()

max = []
for i in range(len(objs)):
    obj = objs[i]
    field = fields[i]
    pre = folder + '/' + obj + '_' + field + '_' + file + '_'
    print(pre)
    base = '_out.pkl'
    extra = pre + 'extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res = pre + 'res' + base
    sed = pre + 'sed' + base
    restwave = pre + 'restwave' + base
    spec = pre + 'spec' + base
    spswave = pre + 'spswave' + base
    chisq = pre + 'chisq' + base
    justchi = pre + 'justchi' + base

    files = [extra, res, sed, restwave, spec, spswave, chisq, justchi]
    print(files)
    # map.all_plots(files, obj, field, file, loc='upper left', sfh=False, curves=False, sep_uvj=False)

    if os.path.exists(files[1]) and os.path.exists(files[7]):
        print('hi')
        if sfh:
            print('sfh')
            print(files[0])
            print_sfh.plotter(files[0], specific=True)
        elif chi_stuff:
            print('chi_stuff')
            with open(files[1], 'rb') as res:
                results = pickle.load(res)
            mask = results['obs']['phot_mask']  # mask out -99.0 values
            with open(files[7], 'rb') as justchi:
                chi = pickle.load(justchi)
            max.append(abs(chi[mask]).max())
            print(1)
        else:
            # PLOT SEPARATE UVJ
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

            with open(files[7], 'rb') as justchi:
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
            chi = chi[mask]

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

            fig = plt.figure()
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
            #    gs = gridspec.GridSpec(1, 2)
            # prepares fig to hold two sets side-by-side of: two axes, one atop other, size ratio 3:1
            ax1 = plt.subplot(gs[0])  # ax1 = bigger upper axis
            # ax3 = plt.subplot(gs[1], sharex=ax1)  # ax3 = smaller lower axis
            ax1.set_yscale("log")
            ax1.set_xscale("log")

            ax1.plot(sps_wave, spec_jan, color='k', alpha=0.5, label=r'Model Spectrum')  # plot spectrum
            ax1.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10,
                     markeredgewidth=1.25,
                     markeredgecolor='k', label=r'Model Photometry')  # plot best fit model
            ax1.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='purple',
                         label=r'Observed Photometry')  # plot observations
            ax1.set_ylabel(r'Flux [$\mu$Jy]', fontsize=fs_text)  # 30

            # ax1.text(1000, 50, 'EELG', fontsize=20)
            # ax1.text(textx, texty, 'COSMOS-1824 (EELG)', fontsize=20)  # 700, 600
            ax1.text(textx, texty, field + '-' + obj, fontsize=20)  # 700, 600
            plt.setp(ax1.get_xticklabels(), visible=False)  # hide xtick labels on upper axis
            # ax1.set_xlabel(r'Rest frame wavelength')
            ax3 = plt.subplot(gs[1], sharex=ax1)  # ax2 = smaller lower axis
            ax3.plot(wave_rest, chi, 'o', color='k')  # plot chi
            plt.axhline(y=0, color='k')  # plot horizontal line at y=0 on chi plot
            yticks = ax3.yaxis.get_major_ticks()  # show ytick labels on lower axis
            yticks[-1].label1.set_visible(False)  # hide uppermost ytick label on lower axis to prevent overlap
            ax3.set_ylabel(r'$\chi$', fontsize=fs_text)
            # ax3.set_xlabel(r'Rest frame wavelength', fontsize=fs)  # 30
            ax1.legend(numpoints=1, loc='upper left', prop={'size': 20})  # , line2) ... , r'$\chi$']
            # plt.subplots_adjust(hspace=.0)

            loc_uvj = 1
            inset_axes(ax1, width=8 * 0.32, height=8 * 0.28, loc=loc_uvj)
            # NOTE: height = 0.875 width (3.5 units vs 4 units) 20%/0.875=0.22857 --> if height 28%, set width to 32%
            # create inset axis: width (%), height (inches), location
            # loc=1 (upper right), loc=2 (upper left) --> loc=3 (lower left), loc=4 (lower right); loc=7 (center right)
            # https://stackoverflow.com/questions/10824156/matplotlib-legend-location-numbers
            uvj.uvj_plot(-1, 'all', objlist=[obj + '_' + field], title=False, labels=False, lims=True,
                         size=20, show=False, col=['purple'])
            print('uvj1')

            plt.subplots_adjust(wspace=.0, hspace=.0)

            xmin = 10 ** 3
            xmax = 2.5 * 10 ** 4
            ymin = 5*10**-3  # 10**-2  # 10 ** -1  # 6*10**-2
            ymax = 5*10**2  # 8*10**3  # 4 * 10 ** 3  # 10**4
            ax1.set_xlim(xmin, xmax)  # 700, xmax
            ax1.set_ylim(ymin, ymax)
            ax3.set_ylim(-2.99, 2.99)
            fs_ticks = 25
            ax1.set_xticks([10 ** 3, 2 * 10 ** 3, 5 * 10 ** 3, 10 ** 4, 2 * 10 ** 4])  # technically works
            ax1.set_xticklabels([r'$10^3$', r'$2\times10^3$', r'$5 \times 10^3$', r'$10^4$', r'$2\times10^4$'],
                                size=fs_ticks)
            ax1.set_yticks([10**-2, 10**-1, 10**0, 10**1, 10**2])  # technically works
            ax1.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$'], size=fs_ticks)
            ax3.set_yticks([-2, -1, 0, 1, 2])  # technically works
            ax3.set_yticklabels([r'$-2$', r'$-1$', r'$0$', r'$1$', r'$2$'], size=fs_ticks)

            ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs)
            ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
            ax3.tick_params('x', length=3, width=1, which='both', labelsize=fs)
            ax3.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)

            # xlabel!
            fig.text(0.5, 0.03, r'Wavelength (Rest) [$\rm \AA$]', ha='center', va='bottom', fontsize=fs_text)  # 30
            plt.show()


if chi_stuff:
    print(np.percentile(max, [16, 50, 84]))
