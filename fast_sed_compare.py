import os
import make_all_plots as map
import pickle
import numpy as np
import print_sfh
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import uvj
import stellar_ages as sa
import get_mass_dust as gmd

home = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
ff = 'fast_fits/'
cdfs_04 = ff + 'cdfs_Z004_EL/'
cdfs_02 = ff + 'cdfs_Z02_EL/'
uds_04 = ff + 'uds_Z004_EL/'
uds_02 = ff + 'uds_Z02_EL/'
cosmos_04 = ff + 'cosmos_Z004_EL/'
cosmos_02 = ff + 'cosmos_Z02_EL/'
cdfs_pre = 'cdfs.v1.6.9_'
cosmos_pre = 'cosmos.v1.3.6_'
uds_pre = 'uds.v1.5.8_'

eelgs = 1  # True=1
folders = ['pkl_efifty/', 'pkl_nvary/']  # ['pkl_efifty', 'pkl_nvary']  # ['pkl_evar', 'pkl_nvary']  # ['pkl_emask', 'pkl_nmask']
folder = 'pkl_efifty/'
out_fold = 'out_efifty/'
file = 'fifty'  # 'fifty'  # 'vary'  # 'newmask'
base = ['fifty', 'vary']
fs = 20  # 30
fs_text = 30
textx = 1100
texty = 150  # 500

objs = []
fields = []
count = 0

eelg_objs, eelg_fields, lbg_objs, lbg_fields = sa.get_gal_lists(base, objlists=True)
full_set = [[eelg_objs, eelg_fields], [lbg_objs, lbg_fields]]

'''
for place in [cdfs_04, cosmos_04, uds_04]:
    for file in os.listdir(place):
        for i in range(len(eelg_objs)):

for i in range(len(eelg_objs)):
    print(eelg_objs[i])
'''

for [objs, fields] in full_set:
    print([objs, fields])
    for i in range(len(objs)):
        obj = objs[i]
        field = fields[i]
        pre = folder + '/' + str(obj) + '_' + field + '_' + file + '_'
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

        if os.path.exists(files[1]):
            print('exist pkl')
            hjon = '/home/jonathan/'
            mass = hjon + 'Comp_10_zm_EL_Z004.dat'  # home + 'Comp_10_zm_EL_Z002.dat'  # home + 'Comp_10_zm_ZFOURGE.dat'
            if fields[i] == 'cosmos':
                place = cosmos_04
                pre = cosmos_pre
                main = '/home/jonathan/cosmos/cosmos.v1.3.8.cat'
                zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'  # redshift catalog
                fc_ob = int(eelg_objs[i]) + 100000
            elif fields[i] == 'cdfs':
                place = cdfs_04
                pre = cdfs_pre
                main = '/home/jonathan/cdfs/cdfs.v1.6.11.cat'
                zname = '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout'
                fc_ob = int(eelg_objs[i])
            elif fields[i] == 'uds':
                place = uds_04
                pre = uds_pre
                main = '/home/jonathan/uds/uds.v.5.10.cat'
                zname = '/home/jonathan/uds/uds.v1.5.8.zout'
                fc_ob = int(eelg_objs[i]) + 200000

            main = np.loadtxt(main)
            zname = np.loadtxt(zname)

            for l in range(len(main)):
                if int(main[l][0]) == eelg_objs[i]:
                    print('yep')
                    if main[l][-1] >= 0:  # if z_spec exists for this galaxy
                        red = main[l][-1]  # index good for all main cats
                    else:
                        red = zname[l][17]  # using zout catalog, z_peak = z_phot; index good for all zout cats

            fast_file = place + pre + str(eelg_objs[i]) + '.awk.fit'
            if os.path.exists(fast_file):

                # GET FAST MASS
                with open(mass) as tm:
                    for line in tm:
                        if line[0] != '#':
                            cols = line.split()
                            print(cols[0], str(fc_ob))
                            if cols[0] == str(fc_ob):
                                fast_mass = cols[2]
                print(fast_mass)
                # GET PROSP MASS
                out = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + out_fold
                for plop in os.listdir(out):
                    if plop.startswith(str(eelg_objs[i])) and plop.endswith(".h5"):
                        get = gmd.printer(out + plop)
                        prosp_mass = get[0]
                print(prosp_mass)
                print('these!')
                print('exist fast')
                print(fast_file)
                with open(fast_file, 'r') as ffile:
                    wls = []
                    flx = []
                    for line in ffile:
                        if line[0] != '#':
                            cols = line.split()
                            wls.append(float(cols[0]) / (1 + red))
                            # 1 Jy = 10-23 erg s-1 cm-2 Hz-1 --> erg * 10^-19
                            # so 1 erg s-1 cm-2 Hz-1 = 10^23 Jy
                            # so 10^-19 erg s-1 cm-2 Hz-1 = 10^4 Jy
                            # ergs
                            flam = float(cols[1]) * 10 ** -19  # erg s^-1 cm^-2 Angstrom^-1
                            # flx in F_lambda --> F_nu = F_lam * lam^2 / c = erg s^-1 cm^-2 Hz-1
                            # erg s-1 cm-2 A-1 * A^2 * A-1 * Hz-1 = erg cm-2 s-1 Hz-1
                            fnu = flam * (float(cols[0]) / (1 + red))**2 / (2.998*10**18)
                            fjy = fnu * 10**23  # 1 Jy = 10^-23 erg s-1 cm-2 Hz-1 --> f[erg,etc.]*10^23 Jy/erg = f[jy]
                            fujy = fjy * 10**6  # Jy to uJy
                            flx.append(fujy)
                print(np.median(flx))

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
                    spec_jan.append(spec[i] * factor)

                fig = plt.figure()
                gs = gridspec.GridSpec(1, 1)  # , height_ratios=[3, 1])
                #    gs = gridspec.GridSpec(1, 2)
                # prepares fig to hold two sets side-by-side of: two axes, one atop other, size ratio 3:1
                ax1 = plt.subplot(gs[0])  # ax1 = bigger upper axis
                # ax3 = plt.subplot(gs[1], sharex=ax1)  # ax3 = smaller lower axis
                ax1.set_yscale("log")
                ax1.set_xscale("log")

                ax1.plot(wls, flx, color='r', alpha=0.5, label=r'FAST model fit')
                ax1.plot(sps_wave, spec_jan, color='k', alpha=0.5, label=r'Model Spectrum')  # plot spectrum
                ax1.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10,
                         markeredgewidth=1.25,
                         markeredgecolor='k', label=r'Model Photometry')  # plot best fit model
                ax1.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='purple',
                             label=r'Observed Photometry')  # plot observations
                ax1.set_ylabel(r'Flux [$\mu$Jy]', fontsize=fs_text)  # 30

                # ax1.text(1000, 50, 'EELG', fontsize=20)
                # ax1.text(textx, texty, 'COSMOS-1824 (EELG)', fontsize=20)  # 700, 600
                ax1.text(textx, texty, field.upper() + '-' + str(obj), fontsize=20)  # 700, 600
                ax1.text(textx, 70, r'Prospector Mass [log$_{10}$(M$_\odot$)]: ' + str(prosp_mass), fontsize=20)
                ax1.text(textx, 30, r'FAST Mass [log$_{10}$(M$_\odot$)]: ' + str(fast_mass), fontsize=20)
                plt.setp(ax1.get_xticklabels(), visible=False)  # hide xtick labels on upper axis
                # ax1.set_xlabel(r'Rest frame wavelength')
                '''
                ax3 = plt.subplot(gs[1], sharex=ax1)  # ax2 = smaller lower axis
                ax3.plot(wave_rest, chi, 'o', color='k')  # plot chi
                plt.axhline(y=0, color='k')  # plot horizontal line at y=0 on chi plot
                yticks = ax3.yaxis.get_major_ticks()  # show ytick labels on lower axis
                yticks[-1].label1.set_visible(False)  # hide uppermost ytick label on lower axis to prevent overlap
                ax3.set_ylabel(r'$\chi$', fontsize=fs_text)
                # ax3.set_xlabel(r'Rest frame wavelength', fontsize=fs)  # 30
                # plt.subplots_adjust(hspace=.0)
                '''
                ax1.legend(numpoints=1, loc='upper left', prop={'size': 20})  # , line2) ... , r'$\chi$']

                loc_uvj = 1
                inset_axes(ax1, width=8 * 0.32, height=8 * 0.28, loc=loc_uvj)
                # NOTE: height = 0.875 width (3.5 units vs 4 units) 20%/0.875=0.22857 --> if height 28%, set width to 32%
                # create inset axis: width (%), height (inches), location
                # loc=1 (upper right), loc=2 (upper left) --> loc=3 (lower left), loc=4 (lower right); loc=7 (center right)
                # https://stackoverflow.com/questions/10824156/matplotlib-legend-location-numbers
                uvj.uvj_plot(-1, 'all', objlist=[str(obj) + '_' + field], title=False, labels=False, lims=True,
                             size=20, show=False, col=['purple'])
                print('uvj1')

                plt.subplots_adjust(wspace=.0, hspace=.0)

                xmin = 10 ** 3
                xmax = 2.5 * 10 ** 4
                ymin = 10**-3  # 10 ** -1  # 6*10**-2
                ymax = 4 * 10 ** 3  # 10**4
                ax1.set_xlim(xmin, xmax)  # 700, xmax
                ax1.set_ylim(ymin, ymax)
                # ax3.set_ylim(-2.99, 2.99)
                ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs)
                ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
                # ax3.tick_params('x', length=3, width=1, which='both', labelsize=fs)
                # ax3.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)

                # xlabel!
                fig.text(0.5, 0.03, r'Wavelength (Rest) [$\rm \AA$]', ha='center', va='bottom', fontsize=fs_text)  # 30
                plt.show()


if chi_stuff:
    print(np.percentile(max, [16, 50, 84]))
